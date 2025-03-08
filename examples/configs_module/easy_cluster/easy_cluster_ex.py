# examples/easy_cluster/easy_cluster_ex.py

from pymmseqs.config import EasyClusterConfig

easy_cluster = EasyClusterConfig(
    fasta_files="input.fasta",
    cluster_prefix="output/clusters",
    tmp_dir="output/tmp",
    v=2
)

easy_cluster.run()
