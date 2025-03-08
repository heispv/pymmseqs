# examples/easy_linclust/easy_linclust_ex.py

from pymmseqs.config import EasyLinClustConfig

config = EasyLinClustConfig(
    fasta_files="input.fasta",
    cluster_prefix="output/clusters",
    tmp_dir="output/tmp",
    v=3,
)
config.run()
