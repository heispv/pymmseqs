# examples/createdb/createdb_ex.py

from pymmseqs.config import (
    CreateDBConfig,
    SearchConfig
)

# Create createdb config objects
query_db = CreateDBConfig(
    fasta_file="query.fasta",
    sequence_db="pydb/q_db",
    v=3
)
target_db = CreateDBConfig(
    fasta_file="target.fasta",
    sequence_db="pydb/t_db",
    v=3
)

# Run the config objects
query_db.run()
target_db.run()

# Create search config object
search_config = SearchConfig(
    query_db="pydb/q_db",
    target_db="pydb/t_db",
    alignment_db="pyout/res",
    tmp_dir="tmp_python",
    v=3
)

# Run the search config object
search_config.run()
