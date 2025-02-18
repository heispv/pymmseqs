# examples/easy_cluster/easy_cluster_ex.py

from pymmseqs.config import EasyClusterConfig

easy_cluster = EasyClusterConfig(
    fasta_files="uniprot.fasta",
    cluster_prefix="new_output/clusters",
    tmp_dir="new_output/tmp",
    v=2
)

easy_cluster.run()
