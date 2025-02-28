# pymmseqs/parsers/createdb_parser.py

class CreateDBParser:
    """
    A class for parsing the output of the CreateDBConfig.
    """
    def __init__(self, sequence_db: str):
        self.sequence_db = sequence_db
    
    def to_path(self) -> str:
        """
        Returns the path to the sequence database.

        Returns:
        --------
        str
        """
        return str(self.sequence_db)
