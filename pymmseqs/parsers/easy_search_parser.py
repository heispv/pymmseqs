# pymmseqs/parsers/easy_search_parser.py

import pandas as pd
import csv
from typing import Generator

from ..config import EasySearchConfig

class EasySearchParser:
    """
    A class for parsing the output of the EasySearchConfig.
    """
    def __init__(self, config: EasySearchConfig):
        if not config.format_mode == 4:
            raise ValueError(f"Using EasySearchParser with format_mode={config.format_mode} is not supported. Please use format_mode=4.")
        
        self.alignment_file = config.alignment_file
    
    def to_pandas(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame containing the alignment data.
        """
        return pd.read_csv(self.alignment_file, sep="\t")
    
    def to_list(self) -> list[dict]:
        """
        Returns a list of dictionaries containing the alignment data.
        """
        return self.to_pandas().to_dict(orient="records")
    
    def to_gen(self) -> Generator[dict, None, None]:
        """
        Returns a generator that yields dictionaries for each row in the alignment file.

        Each dictionary represents a row in the TSV file, with keys corresponding to 
        the column names in the header.
        """
        with open(self.alignment_file, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            yield from reader
    
    def to_path(self) -> str:
        """
        Returns a list of file paths for the output files.

        Returns:
        --------
        list of str
        """
        return self.alignment_file
