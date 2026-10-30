"""
The properties of the tsls HMM that sc2ts relies on.

Most of these tests call tsls directly, with the parameters that sc2ts
computes, so that a change in tsls that would break sc2ts fails here with a
clear message. Some also check the same through sc2ts's match_tsls.
"""

import concurrent.futures as cf

import numpy as np
import pytest
import tskit
import tsls
import util

import sc2ts
from sc2ts import inference as si
from sc2ts import jit

NUM_SITES = 40
# The first PANEL_LEFT sites are where node A's mutations are, and the rest
# are where node B's are.
PANEL_LEFT = 10


def reference_ts():
    """The initial sc2ts tree sequence, with the first NUM_SITES sites."""
    tables = util.initial_ts().dump_tables()
    tables.sites.truncate(NUM_SITES)
    tables.nodes.time += 10
    return tables.tree_sequence()


def reference_haplotype(ts):
    alignment = util.reference_array()
    alignment[0] = "A"
    return jit.encode_alleles(alignment)[ts.sites_position.astype(int)]


def other_allele(state):
    return "A" if state != "A" else "C"


def panel_ts(num_b_mutations):
    """
    The reference (node 1) with two children: A (node 2), with a mutation at
    each of the first PANEL_LEFT sites, and B (node 3), with mutations at
    num_b_mutations sites after them.
    """
    ts = reference_ts()
    tables = ts.dump_tables()
    a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=1, metadata={})
    b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=1, metadata={})
    for child in [a, b]:
        tables.edges.add_row(0, ts.sequence_length, parent=1, child=child)
    for site in range(PANEL_LEFT):
        state = other_allele(ts.site(site).ancestral_state)
        tables.mutations.add_row(site=site, node=a, derived_state=state)
    for site in range(PANEL_LEFT, PANEL_LEFT + num_b_mutations):
        state = other_allele(ts.site(site).ancestral_state)
        tables.mutations.add_row(site=site, node=b, derived_state=state)
    tables.sort()
    tables.build_index()
    return tables.tree_sequence()


def node_haplotype(ts, u):
    G = ts.genotype_matrix(
        samples=[u], alleles=tuple(sc2ts.IUPAC_ALLELES), isolated_as_missing=False
    )
    return G[:, 0].astype(np.int8)


def hmm_parameters(ts, num_mismatches, deletions_as_missing=False):
    """
    The (recombination, mismatch) probabilities for tsls, as match_tsls
    computes them.
    """
    num_alleles = 4 if deletions_as_missing else 5
    mu, rho = si.solve_num_mismatches(num_mismatches, num_alleles)
    recombination = tsls.model.recombination_probability(rho, ts.num_nodes)
    odds = tsls.model.mismatch_odds(mu, num_alleles)
    mismatch = tsls.model.mismatch_probability(odds, 5)
    return float(recombination), float(mismatch)


def find_path(ts, h, num_mismatches=3, deletions_as_missing=False, prune_threshold=None):
    """
    Run tsls as match_tsls does, returning (left, right, parent, match), with
    the path in right-to-left order.
    """
    hmm_ts = si.make_hmm_tree_sequence(ts)
    recombination, mismatch = hmm_parameters(ts, num_mismatches, deletions_as_missing)
    m = ts.num_sites
    hmm = tsls.Hmm(
        hmm_ts,
        np.full(m, recombination),
        np.full(m, mismatch),
        prune_threshold=prune_threshold,
    )
    match = np.zeros(m, dtype=np.int8)
    h = np.asarray(h, dtype=np.int8)
    left, right, parent = hmm.find_path(h, 0, m, match)
    return left, right, parent, match


def path_cost(h, parent, match, num_mismatches):
    observed = h != -1
    return num_mismatches * (len(parent) - 1) + int(
        np.sum(h[observed] != match[observed])
    )


class TestAlleles:
    def test_hmm_alleles_are_iupac_codes(self):
        assert si.HMM_ALLELES == list(sc2ts.IUPAC_ALLELES[:5])
        hmm_ts = si.make_hmm_tree_sequence(reference_ts())
        assert np.all(hmm_ts.num_alleles == 5)

    def test_copy_node_exactly(self):
        ts = panel_ts(5)
        for u in [1, 2, 3]:
            h = node_haplotype(ts, u)
            left, right, parent, match = find_path(ts, h)
            _, _, parent = si.copy_node1_for_node0(left, right, parent)
            assert list(parent) == [u]
            np.testing.assert_array_equal(match, h)

    def test_deletion_is_an_allele(self):
        ts = reference_ts()
        tables = ts.dump_tables()
        u = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=1, metadata={})
        tables.edges.add_row(0, ts.sequence_length, parent=1, child=u)
        tables.mutations.add_row(site=3, node=u, derived_state="-")
        tables.sort()
        ts = tables.tree_sequence()
        h = reference_haplotype(ts)
        h[3] = si.DELETION
        _, _, parent, match = find_path(ts, h)
        assert list(parent) == [u]
        assert match[3] == si.DELETION

    def test_missing_data(self):
        ts = panel_ts(5)
        h = node_haplotype(ts, 3)
        h[[0, PANEL_LEFT + 1, NUM_SITES - 1]] = -1
        _, _, parent, match = find_path(ts, h)
        assert list(parent) == [3]
        np.testing.assert_array_equal(match, node_haplotype(ts, 3))

    def test_state_outside_hmm_alleles(self):
        ts = reference_ts()
        tables = ts.dump_tables()
        tables.mutations.add_row(site=3, node=1, derived_state="N")
        with pytest.raises(ValueError, match="not in its site's list of alleles"):
            tsls.TreeSequence(tables.tree_sequence(), alleles=si.HMM_ALLELES)


class TestVestigialRoot:
    @pytest.mark.parametrize("ts", [reference_ts(), panel_ts(5)])
    def test_node0_ties_node1(self, ts):
        # Node 0 has the same Viterbi value as node 1 at every site, so
        # copying node 1 where tsls copies node 0 is exact.
        hmm_ts = si.make_hmm_tree_sequence(ts)
        recombination, mismatch = hmm_parameters(ts, 3)
        m = ts.num_sites
        hmm = tsls.Hmm(hmm_ts, np.full(m, recombination), np.full(m, mismatch))
        h = node_haplotype(ts, 1)
        h[5] = si.DELETION
        _, _, matrix = hmm.viterbi_matrix(h)
        V = matrix.decode([0, 1])
        np.testing.assert_array_equal(V[:, 0], V[:, 1])

    def test_initial_ts(self):
        # Against the initial tree sequence, only nodes 0 and 1 can be copied.
        ts = reference_ts()
        h = reference_haplotype(ts)
        h[2] = si.DELETION
        _, _, parent, _ = find_path(ts, h)
        assert set(parent) <= {0, 1}
        left, right, parent = si.copy_node1_for_node0(
            np.array([30.0, 10.0, 0.0]),
            np.array([40.0, 30.0, 10.0]),
            np.array([1, 0, 1], dtype=np.int32),
        )
        assert list(parent) == [1]
        assert list(left) == [0]
        assert list(right) == [40]

    def test_mutation_on_reference(self):
        tables = reference_ts().dump_tables()
        tables.mutations.add_row(site=0, node=1, derived_state="T")
        with pytest.raises(ValueError, match="vestigial root"):
            si.make_hmm_tree_sequence(tables.tree_sequence())

    def test_other_child_of_node0(self):
        ts = reference_ts()
        tables = ts.dump_tables()
        u = tables.nodes.add_row(time=1, metadata={})
        tables.edges.add_row(0, ts.sequence_length, parent=0, child=u)
        tables.sort()
        with pytest.raises(ValueError, match="vestigial root"):
            si.check_vestigial_root(tables.tree_sequence())


class TestIsolatedNodes:
    def test_isolated_node_not_copied(self):
        # Future nodes are detached before matching, and isolated nodes are
        # not states of the HMM.
        ts = panel_ts(5)
        tables = ts.dump_tables()
        tables.edges.keep_rows(tables.edges.child != 3)
        tables.mutations.keep_rows(tables.mutations.node != 3)
        detached = tables.tree_sequence()
        h = node_haplotype(ts, 3)
        _, _, parent, _ = find_path(detached, h)
        assert 3 not in parent

    def test_future_node_not_copied(self):
        ts = panel_ts(5)
        tables = ts.dump_tables()
        tables.nodes[3] = tables.nodes[3].replace(time=-1)
        tables.sort()
        ts = tables.tree_sequence()
        h = node_haplotype(ts, 3)
        sample = si.Sample("test", "2020-01-01", haplotype=h)
        si.match_tsls([sample], ts, num_mismatches=3, mismatch_threshold=20)
        assert 3 not in sample.hmm_match.parents


class TestCostModel:
    """
    One recombination costs num_mismatches = k mismatches: the query is A on
    the left and B on the right, so it can copy A with j mismatches or switch
    from A to B once. sc2ts's parameters make the switch slightly more
    likely than k mismatches.
    """

    def query(self, ts):
        h = node_haplotype(ts, 2)
        h[PANEL_LEFT:] = node_haplotype(ts, 3)[PANEL_LEFT:]
        return h

    @pytest.mark.parametrize("deletions_as_missing", [False, True])
    @pytest.mark.parametrize("k", [2, 3, 4])
    @pytest.mark.parametrize("delta", [-1, 0, 1])
    def test_switch_or_mismatch(self, k, delta, deletions_as_missing):
        j = k + delta
        ts = panel_ts(j)
        h = self.query(ts)
        left, _, parent, match = find_path(
            ts, h, num_mismatches=k, deletions_as_missing=deletions_as_missing
        )
        if j >= k:
            assert list(parent) == [3, 2]
            assert PANEL_LEFT <= np.searchsorted(ts.sites_position, left[0])
            np.testing.assert_array_equal(match, h)
        else:
            assert list(parent) == [2]
            assert np.sum(match != h) == j

    @pytest.mark.parametrize("mirror", [False, True])
    @pytest.mark.parametrize("delta", [-1, 0, 1])
    def test_match_tsls(self, delta, mirror):
        k = 3
        j = k + delta
        ts = panel_ts(j)
        sample = si.Sample("test", "2020-01-01", haplotype=self.query(ts))
        si.match_tsls(
            [sample],
            ts,
            num_mismatches=k,
            mismatch_threshold=20,
            mirror_coordinates=mirror,
        )
        match = sample.hmm_match
        if j >= k:
            assert match.parents == [2, 3]
            assert match.cost == k
            assert len(match.mutations) == 0
        else:
            assert match.parents == [2]
            assert match.cost == j
            assert len(match.mutations) == j


class TestCostThreshold:
    """
    The prune thresholds that match_tsls uses guarantee the optimal path
    whenever it costs at most the threshold's cost.
    """

    def panel(self, seed):
        rng = np.random.default_rng(seed)
        ts = reference_ts()
        tables = ts.dump_tables()
        parent = 1
        for j in range(8):
            time = 1 - j / 10
            u = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=time, metadata={})
            tables.edges.add_row(0, ts.sequence_length, parent=parent, child=u)
            for site in rng.choice(NUM_SITES, size=4, replace=False):
                state = other_allele(ts.site(site).ancestral_state)
                tables.mutations.add_row(site=site, node=u, derived_state=state)
            parent = rng.choice([1, u])
        tables.sort()
        tables.build_index()
        tables.compute_mutation_parents()
        return tables.tree_sequence()

    @pytest.mark.parametrize("seed", range(5))
    @pytest.mark.parametrize("max_cost", [1, 2, 4])
    def test_optimal_within_cost(self, seed, max_cost):
        k = 3
        ts = self.panel(seed)
        recombination, mismatch = hmm_parameters(ts, k)
        threshold = si.cost_threshold(max_cost, k, mismatch, recombination, ts.num_nodes)
        rng = np.random.default_rng(seed)
        nodes = np.arange(1, ts.num_nodes)
        for _ in range(10):
            a, b = rng.choice(nodes, size=2, replace=False)
            h = node_haplotype(ts, a)
            h[NUM_SITES // 2 :] = node_haplotype(ts, b)[NUM_SITES // 2 :]
            sites = rng.choice(NUM_SITES, size=rng.integers(0, 3), replace=False)
            h[sites] = (h[sites] + 1) % 5
            _, _, parent, match = find_path(ts, h, k)
            optimum = path_cost(h, parent, match, k)
            _, _, parent, match = find_path(ts, h, k, prune_threshold=threshold)
            cost = path_cost(h, parent, match, k)
            assert cost >= optimum
            if optimum <= max_cost:
                assert cost == optimum

    def test_threshold_values(self):
        # With no switches within the cost, the threshold is the odds of
        # max_cost mismatches; a switch costing k mismatches is cheaper than
        # k + 1 mismatches.
        ts = panel_ts(3)
        k = 3
        recombination, mismatch = hmm_parameters(ts, k)
        e = tsls.model.mismatch_odds(mismatch, 5)
        t = si.cost_threshold(2, k, mismatch, recombination, ts.num_nodes)
        assert t == pytest.approx(e**2)
        t = si.cost_threshold(4, k, mismatch, recombination, ts.num_nodes)
        assert t == pytest.approx(e**4)


class TestPathFormat:
    def test_right_to_left_on_site_positions(self):
        ts = panel_ts(4)
        h = node_haplotype(ts, 2)
        h[PANEL_LEFT:] = node_haplotype(ts, 3)[PANEL_LEFT:]
        left, right, parent, match = find_path(ts, h)
        assert list(parent) == [3, 2]
        assert left[-1] == 0
        assert right[0] == ts.sequence_length
        np.testing.assert_array_equal(left[:-1], right[1:])
        assert left[0] in ts.sites_position
        # The switch is between the last site that needs A and the first
        # that needs B.
        x = np.searchsorted(ts.sites_position, left[0])
        assert PANEL_LEFT <= x <= PANEL_LEFT + 4 - 1 + 1

    def test_matched_haplotype_whole_sequence(self):
        ts = panel_ts(4)
        h = node_haplotype(ts, 3)
        _, _, _, match = find_path(ts, h)
        assert np.all(match >= 0)


class TestThreads:
    def test_shared_tree_sequence(self):
        ts = panel_ts(5)
        hmm_ts = si.make_hmm_tree_sequence(ts)
        recombination, mismatch = hmm_parameters(ts, 3)
        m = ts.num_sites
        rng = np.random.default_rng(1)
        H = rng.integers(0, 5, size=(40, m)).astype(np.int8)

        def run(h):
            hmm = tsls.Hmm(hmm_ts, np.full(m, recombination), np.full(m, mismatch))
            match = np.zeros(m, dtype=np.int8)
            path = hmm.find_path(h, 0, m, match)
            assert hmm.total_memory > 0
            return [array.tolist() for array in path] + [match.tolist()]

        serial = [run(h) for h in H]
        with cf.ThreadPoolExecutor(max_workers=4) as executor:
            threaded = list(executor.map(run, H))
        assert serial == threaded
        assert serial == [run(h) for h in H]
