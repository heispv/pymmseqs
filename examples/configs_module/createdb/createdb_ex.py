# examples/createdb/createdb_ex.py

from pymmseqs.config import CreateDBConfig

# Create and validate config
config = CreateDBConfig(
    fasta_file="input.fasta",
    sequence_db="output/db"
)

# To use a config file, uncomment the following line
# config = CreateDBConfig.from_yaml("config.yaml")

config.run()
