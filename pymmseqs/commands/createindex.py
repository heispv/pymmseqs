# pymmseqs/commands/createindex.py

from pathlib import Path
from typing import Union

from ..config import CreateIndexConfig
from ..parsers import CreateIndexParser
from ..utils import tmp_dir_handler

def createindex(
    sequence_db: Union[str, Path],
    tmp_dir: Union[str, Path] = None,

    s: float = 7.5,
    k: int = 0,
    v: int = 3,
    threads: int = 14,
    compressed: bool = False,
    create_lookup: int = 0,
    
    search_type: int = 0,
    headers_split_mode: int = 0,
    max_seqs: int = 300,
    max_seq_len: int = 65535,
) -> CreateIndexParser:
    """
    Create an index for a sequence database using MMseqs2.
    
    Parameters
    ----------
    `sequence_db` : Union[str, Path]
        Path to MMseqs2 sequence database created with createdb.
        
    `tmp_dir` : Union[str, Path]
        Temporary directory for intermediate files.
        If not provided, a temporary directory will be created in the same directory as the sequence_db.
    
    `s` : float, optional
        Sensitivity.
        - Options: 1.0 (faster), 4.0 (fast), 7.5 (sensitive)
        - Default: 7.5
        
    `k` : int, optional
        k-mer length.
        - 0: automatically set to optimum (default)
        
    `v` : int, optional
        Verbosity level (0-3).
        - Default: 3
        
    `threads` : int, optional
        Number of threads to use.
        - Default: 14
        
    `compressed` : bool, optional
        Use compressed database.
        - True
        - False (default)
        
    `create_lookup` : int, optional
        Create lookup file.
        - 0: no lookup file (default)
        - 1: create lookup file
        
    `search_type` : int, optional
        Search type
        - 0: auto (default)
        - 1: amino acid
        - 2: nucleotide
        - 3: translated
        
    `headers_split_mode` : int, optional
        Header split mode
        - 0: split position (default)
        - 1: original header
        
    `max_seqs` : int, optional
        Maximum results per query sequence passing prefilter.
        - Default: 300
        
    `max_seq_len` : int, optional
        Maximum sequence length.
        - Default: 65535
    
    Returns
    -------
    CreateIndexParser
        Parser for the created index.
    """

    tmp_dir = tmp_dir_handler(
        tmp_dir=tmp_dir,
        output_file_path=sequence_db
    )

    config = CreateIndexConfig(
        sequence_db=sequence_db,
        tmp_dir=tmp_dir,
        s=s,
        k=k,
        v=v,
        threads=threads,
        compressed=compressed,
        create_lookup=create_lookup,
        search_type=search_type,
        headers_split_mode=headers_split_mode,
        max_seqs=max_seqs,
        max_seq_len=max_seq_len,
    )

    config.run()

    return CreateIndexParser(config)
