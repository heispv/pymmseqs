# examples/createtaxdb/createtaxdb_ex.py

from pymseqs.commands import createdb, createtaxdb

createdb("input.fasta", "output/db", tmp_dir="output/tmp")

createtaxdb(
    sequence_db="output/db",
    tmp_dir="output/tmp"
)
