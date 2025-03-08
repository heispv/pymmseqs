# examples/cluster/cluster_ex.py

from pymseqs.commands import createdb
from pymseqs.commands import cluster

createdb("input.fasta", "output/db_output/db")

cluster(
    "output/db_output/db",
    "output/cluster_output/clusters",
    "output/tmp"
)

# TODO: Then you should also use createtsv to create a TSV file from the clusters (human readable format)
# createtsv("output/db", "output/clusters", "output/clusters.tsv")
