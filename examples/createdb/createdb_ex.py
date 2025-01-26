# examples/createdb/createdb_ex.py

from pymseqs.commands import createdb

createdb("input.fasta", "output/db", write_lookup=False, v=3)
