# examples/easy_cluster/easy_cluster_ex.py

from pymseqs.commands import easy_cluster

easy_cluster(
    "input.fasta",
    "output/clusters",
    "output/tmp",
    v=3,
)
