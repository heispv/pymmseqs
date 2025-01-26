# examples/createdb/createdb_ex.py

from pymseqs.commands import createdb, search

createdb("target.fasta", "pydb/t_db", v=3)
createdb("query.fasta", "pydb/q_db", v=3)
search ("pydb/q_db", "pydb/t_db", "pyout/res", tmp_dir="tmp_python", v=3)
