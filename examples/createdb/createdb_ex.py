# examples/createdb/createdb_ex.py

from pymseqs.commands import createdb
import os


createdb('input.fasta', 'db', write_lookup=0, verbosity=3)