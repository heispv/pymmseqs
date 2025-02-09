# examples/easy_taxonomy/easy_taxonomy_ex.py

from pymseqs.commands import createdb, createtaxdb, easy_taxonomy

createdb("target.fasta", "pydb/t_db", tmp_dir="tmp")
createdb("query.fasta", "pydb/q_db", tmp_dir="tmp")
createtaxdb("pydb/t_db", tmp_dir="tmp")
easy_taxonomy("pydb/q_db", "pydb/t_db", "pyout/taxonomy_res", tmp_dir="tmp")
