# pymmseqs/parsers/createindex_parser.py

from ..config import CreateIndexConfig

class CreateIndexParser:
    """
    A class for parsing the output of the CreateIndexConfig.
    """
    def __init__(self, config: CreateIndexConfig):
        self.sequence_db = config.sequence_db
    
    def to_path(self) -> str:
        """
        Returns the path to the sequence database.

        Returns:
        --------
        str
        """
        return str(self.sequence_db)
