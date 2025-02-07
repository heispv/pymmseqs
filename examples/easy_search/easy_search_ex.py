# examples/easy_search/easy_search_ex.py

from pymseqs.commands import easy_search

easy_search (
    query_fasta="query.fasta",
    target_fasta_or_db="target.fasta",
    alignment_file="output/result",
    tmp_dir="output/tmp_python",
    format_mode=4
    )
