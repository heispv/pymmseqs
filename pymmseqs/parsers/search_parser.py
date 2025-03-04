# pymmseqs/parsers/search_parser.py

import pandas as pd
import csv
from typing import Generator

from ..config import (
    SearchConfig,
    ConvertAlisConfig
)
class SearchParser:
    """
    A class for parsing the output of the SearchConfig.
    """
    def __init__(self, config: SearchConfig):
        self.query_db = config.query_db
        self.target_db = config.target_db
        self.alignment_db = config.alignment_db

        self._readable = False
    
    def _run_convertalis(self) -> None:
        """
        Runs the convertalis command to convert the alignment database to a readable format.
        """
        print("Output is not readable. Executing convertalis command to convert the alignment database to a readable format.")

        config = ConvertAlisConfig(
            query_db=self.query_db,
            target_db=self.target_db,
            alignment_db=self.alignment_db,
            alignment_file=f"{self.alignment_db}.tsv",
            format_mode=4
        )
        config.run()
        self._readable = True
        
    def to_pandas(self) -> pd.DataFrame:
        """
        Returns a pandas DataFrame containing the alignment data.
        """
        if not self._readable:
            self._run_convertalis()
        
        return pd.read_csv(f"{self.alignment_db}.tsv", sep="\t")
    
    def to_list(self) -> list[dict]:
        """
        Returns a list of dictionaries containing the alignment data.
        """
        if not self._readable:
            self._run_convertalis()
        
        return self.to_pandas().to_dict(orient="records")
    
    def to_gen(self) -> Generator[dict, None, None]:
        """
        Returns a generator that yields dictionaries for each row in the alignment file.

        Each dictionary represents a row in the TSV file, with keys corresponding to 
        the column names in the header.
        """
        if not self._readable:
            self._run_convertalis()
        
        with open(f"{self.alignment_db}.tsv", 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            yield from reader
    
    def to_path(self) -> str:
        """
        Returns a list of file paths for the output files.

        Returns:
        --------
        list of str
        """
        
        return f"{self.alignment_db}"
