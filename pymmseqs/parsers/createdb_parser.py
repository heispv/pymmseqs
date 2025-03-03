# pymmseqs/parsers/createdb_parser.py

from ..config import CreateDBConfig

class CreateDBParser:
    """
    A class for parsing the output of the CreateDBConfig.
    """
    def __init__(self, config: CreateDBConfig):
        self.sequence_db = config.sequence_db
    
    def to_path(self) -> str:
        """
        Returns the path to the sequence database.

        Returns:
        --------
        str
        """
        return str(self.sequence_db)
