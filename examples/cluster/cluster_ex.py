# examples/cluster/cluster_ex.py

from pymseqs.commands import createdb
from pymseqs.commands import cluster

createdb("input.fasta", "output/db")

cluster(
    "output/db",
    "output/clusters",
    "output/tmp",
    min_seq_id=1)

# TODO: Then you should also use createtsv to create a TSV file from the clusters (human readable format)
# createtsv("output/db", "output/clusters", "output/clusters.tsv")
