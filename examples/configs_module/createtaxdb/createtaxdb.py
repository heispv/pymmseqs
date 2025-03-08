# examples/createtaxdb/createtaxdb_ex.py

from pymmseqs.config import (
    CreateDBConfig,
    CreateTaxDBConfig
)

create_db = CreateDBConfig(
    sequence_db="input.fasta",
    db_name="output/db",
    tmp_dir="output/tmp"
)
create_db.run()

create_tax_db = CreateTaxDBConfig(
    sequence_db="output/db",
    tmp_dir="output/tmp"
)
create_tax_db.run()
