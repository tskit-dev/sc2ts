import numpy as np
import tskit

import sc2ts
from sc2ts import data_import
from sc2ts import inference as si


def reference_sequence():
    """The built-in SARS-CoV-2 reference, with "X" prepended (1-based)."""
    _, sequence = data_import.read_fasta(data_import.data_path / "reference.fasta")
    return sequence


def reference_array():
    """The built-in SARS-CoV-2 reference as a char array with "X" at index 0."""
    return np.array(list(reference_sequence()))


def initial_ts(problematic_sites=None):
    """
    Build an initial tree sequence using the built-in SARS-CoV-2 reference.
    This mirrors what the CLI does by default (see ``sc2ts.cli.get_reference``).
    """
    return si.initial_ts(
        reference_sequence=reference_sequence(),
        reference_id=sc2ts.REFERENCE_STRAIN,
        reference_date=sc2ts.REFERENCE_DATE,
        problematic_sites=problematic_sites,
    )


def get_match_db(ts, db_path, samples, date, num_mismatches):
    si.MatchDb.initialise(db_path)
    match_db = si.MatchDb(db_path)
    match_db.add(samples, date, num_mismatches)
    match_db.create_mask_table(ts)
    return match_db


def example_binary(n, date="2020-01-01"):
    base = initial_ts()
    tables = base.dump_tables()
    tree = tskit.Tree.generate_balanced(n, span=base.sequence_length)
    binary_tables = tree.tree_sequence.dump_tables()
    binary_tables.nodes.time += 1
    tables.nodes.time += np.max(binary_tables.nodes.time) + 1
    binary_tables.edges.child += len(tables.nodes)
    binary_tables.edges.parent += len(tables.nodes)
    for j, node in enumerate(binary_tables.nodes):
        md = {}
        if node.flags == tskit.NODE_IS_SAMPLE:
            md["strain"] = f"x{j}"
            md["date"] = date
        tables.nodes.append(node.replace(metadata=md))
    for edge in binary_tables.edges:
        tables.edges.append(edge)
    # FIXME brittle
    tables.edges.add_row(0, base.sequence_length, parent=1, child=tree.root + 2)
    tables.sort()
    return tables.tree_sequence()
