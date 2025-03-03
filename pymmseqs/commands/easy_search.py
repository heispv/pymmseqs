# pymmseqs/commands/easy_search.py

from pathlib import Path
from typing import Union, List

from ..config import EasySearchConfig
from ..parsers import EasySearchParser

def easy_search(
    query_fasta: Union[str, Path, List[Union[str, Path]]],
    target_fasta_or_db: Union[str, Path],
    alignment_file: Union[str, Path],
    tmp_dir: Union[str, Path],

    # Optional parameters
    s: float = 5.7,
    e: float = 0.001,
    min_seq_id: float = 0.0,
    c: float = 0.0,
    max_seqs: int = 300,

) -> EasySearchParser:
    """
    Required parameters
    ----------
    `query_fasta` : Union[str, Path]
        Path to one or more query FASTA files. Can be compressed with .gz or .bz2.

    `target_fasta_or_db` : Union[str, Path]
        Path to a target FASTA file (optionally compressed) or an MMseqs2 target database.

    `alignment_file` : Union[str, Path]
        Path to the output file where alignments will be stored.

    `tmp_dir` : Union[str, Path]
        Temporary directory for intermediate files. Will be created if not existing.
    
    Optional parameters
    -------------------
    `s` : float, optional
        Sensitivity
        - 1.0: faster
        - 4.0: fast
        - 5.7 (default)
        - 7.5: sensitive
    
    `e` : float, optional
        E-value threshold (range 0.0, inf)
        - 0.001 (default)
    
    `min_seq_id` : float, optional
        Minimum sequence identity (range 0.0, 1.0)
        - 0.0 (default)
    
    `c` : float, optional
        Coverage threshold for alignments
        - 0.0 (default)
        - Determines the minimum fraction of aligned residues required for a match, based on the selected cov_mode
    
    `max_seqs` : int, optional
        Maximum results per query passing prefilter
        - 300 (default)
        - Higher values increase sensitivity but may slow down the search
    
    Returns
    -------
    EasySearchParser object
        - An EasySearchParser instance that provides methods to access and parse the alignment data.
    """

    config = EasySearchConfig(
        query_fasta=query_fasta,
        target_fasta_or_db=target_fasta_or_db,
        alignment_file=alignment_file,
        tmp_dir=tmp_dir,
        s=s,
        e=e,
        min_seq_id=min_seq_id,
        c=c,
        max_seqs=max_seqs,
        format_mode=4
    )

    config.run()

    return EasySearchParser(config)
