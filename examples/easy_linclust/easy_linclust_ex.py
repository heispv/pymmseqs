# examples/easy_linclust/easy_linclust_ex.py

from pymseqs.commands import easy_linclust

easy_linclust(
    "input.fasta",
    "output/clusters",
    "output/tmp",
    v=3,
)
