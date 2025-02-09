# examples/easy_taxonomy/easy_taxonomy_ex.py

from pymseqs.commands import createdb, easy_taxonomy

createdb("target.fasta", "pydb/t_db")
createdb("query.fasta", "pydb/q_db")
# Right now we don't have the createtaxdb command, so it doens't work
createtaxdb("pydb/t_db", tmp_dir="tmp")
easy_taxonomy("pydb/q_db", "pydb/t_db", "pyout/taxonomy_res", tmp_dir="tmp")
