# pymmseqs/commands/easy_search.py

from pathlib import Path
from typing import Union, List
import pandas as pd
from ..config import EasySearchConfig

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

    # Output parameters
    output_mode: str = None

) -> Union[pd.DataFrame, List[dict], Path, None]:
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
    
    Output parameters
    ------------------
    `output_mode` : str, optional
        Output mode
        - "to_pandas": return a pandas DataFrame
        - "to_dict": return a list of dictionaries
        - "to_rel_path": return a relative path to the alignment file
        - "to_abs_path": return an absolute path to the alignment file
    
    Returns
    -------
    pd.DataFrame, List[dict], Path, or None
        - Based on the `output_mode` parameter
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

    if output_mode == "to_pandas":
        return pd.read_csv(alignment_file, sep="\t")
    elif output_mode == "to_dict":
        return pd.read_csv(alignment_file, sep="\t").to_dict(orient="records")
    elif output_mode == "to_rel_path":
        return str(alignment_file)
    elif output_mode == "to_abs_path":
        return config.alignment_file
    else:
        return None
