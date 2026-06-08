# sc2ts <img align="right" width="145" height="90" src="https://raw.githubusercontent.com/tskit-dev/sc2ts/main/docs/sc2ts.png">

Sc2ts stands for "SARS-CoV-2 tree sequence" (pronounced "scoots" optionally) and consists of

1. A method to infer Ancestral Recombination Graphs (ARGs) from SARS-CoV-2 genome sequence data at pandemic scale.
2. A lightweight wrapper around [tskit Python APIs](https://tskit.dev/tskit/docs/stable/python-api.html) specialised for the output of sc2ts which enables efficient node metadata access.
3. A lightweight wrapper around [Zarr Python](https://zarr.dev) which enables convenient and efficient access to the full Viridian dataset (alignments and metadata) in a single file using the [VCF Zarr specification](https://doi.org/10.1093/gigascience/giaf049).

For details on the software, please see the online [documentation](https://tskit.dev/sc2ts/docs). For information on the method and an inferred ARG, please see this [preprint](https://www.biorxiv.org/content/10.1101/2023.06.08.544212v3):

> Shing H. Zhan, Yan Wong, Anastasia Ignatieva, Katherine Eaton, Isobel Guthrie, Benjamin Jeffery, Duncan S. Palmer, Carmen Lia Murall, Sarah P. Otto, and Jerome Kelleher (2025) _A Pandemic-Scale Ancestral Recombination Graph for SARS-CoV-2_. bioRxiv: 2023.06.08.544212; doi: https://doi.org/10.1101/2023.06.08.544212
