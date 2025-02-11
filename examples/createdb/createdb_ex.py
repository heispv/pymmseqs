# examples/createdb/createdb_ex.py

from pymseqs.config import CreateDBConfig
from pymseqs.commands import createdb

# Create and validate config
config = CreateDBConfig(
    input_fasta="input.fasta",
    db_name="output/db"
)

# To use a config file, uncomment the following line

# config = CreateDBConfig.from_yaml("config.yaml")

createdb(config)
