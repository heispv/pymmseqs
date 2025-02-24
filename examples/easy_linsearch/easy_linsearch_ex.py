# examples/easy_linsearch/easy_linsearch_ex.py

from pymmseqs.config import EasyLinSearchConfig

config = EasyLinSearchConfig(
    query_fasta="query.fasta",
    target_fasta_or_db="target.fasta",
    alignment_file="output/result",
    tmp_dir="output/tmp_python"
    )

config.run()
