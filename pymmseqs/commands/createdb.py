# pymmseqs/commands/createdb.py

from pathlib import Path
from typing import Union, List

from ..config import CreateDBConfig
from ..parsers import CreateDBParser

def createdb(
    fasta_file: Union[List[Union[str, Path]], Union[str, Path]],
    sequence_db: Union[str, Path],

    # Optional parameters
    shuffle: bool = True,
    compressed: bool = False,
    createdb_mode: int = 0,
    dbtype: int = 0,

) -> CreateDBParser:
    """
    Create a MMseqs2 database from a FASTA file and save it to the specified path prefix

    Parameters
    ----------
    `fasta_file` : Union[List[Union[str, Path]], Union[str, Path]]
    Path(s) to the input FASTA file(s). This can be:
        - A single string or Path object (e.g., `"input.fasta"` or `Path("input.fasta")`)
        - A list of strings or Path objects (e.g., `["input1.fasta", "input2.fasta"]` or 
          `[Path("input1.fasta"), Path("input2.fasta")]`)
    
    `sequence_db` : Union[str, Path]
        Database path prefix, including the desired directory structure (e.g., `"output/dbs/mydb"`)
    
    `shuffle` : bool, optional
        Shuffle the input database entries
        - True (default)
        - False
    
    `compressed` : bool, optional
        Compress the output files
        - True
        - False (default)
    
    `createdb_mode` : int, optional
        Database creation mode
        - 0: Copy data (default)
        - 1: Soft-link data and write a new index (only works with single-line FASTA/Q)
    
    `dbtype` : int, optional
        Database type
        - 0: Auto-detect (default)
        - 1: Amino acid sequences
        - 2: Nucleotide sequences
    
    Returns
    -------
    CreateDBParser object
        - An CreateDBParser instance that provides methods to access and parse the sequence database.
    """

    config = CreateDBConfig(
        fasta_file=fasta_file,
        sequence_db=sequence_db,
        shuffle=shuffle,
        compressed=compressed,
        createdb_mode=createdb_mode,
        dbtype=dbtype,
    )

    config.run()

    return CreateDBParser(config)
