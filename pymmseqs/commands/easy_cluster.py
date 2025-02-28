# pymmseqs/commands/easy_cluster.py

import os
from pathlib import Path
from typing import Union, List, Generator

from ..config import EasyClusterConfig
from ..parsers import EasyClusterParser

def easy_cluster(
    # Required parameters
    fasta_files: Union[str, Path, List[Union[str, Path]]],
    cluster_prefix: Union[str, Path],
    tmp_dir: Union[str, Path],

    # Optional parameters
    min_seq_id: float = 0.0,
    s: float = 4.0,
    c: float = 0.8,
    cov_mode: int = 0,
    e: float = 0.001,
    cluster_mode: int = 0,

) -> 'EasyClusterParser':
    """
    Perform sequence clustering using MMseqs2.

    Parameters
    ----------
    `fasta_files` : Union[str, Path, List[Union[str, Path]]]
        Path to MMseqs2 sequence database created with createdb.

    `cluster_prefix` : Union[str, Path]
        Output cluster database path prefix (will create multiple files with this prefix).

    `tmp_dir` : Union[str, Path]
        Temporary directory for intermediate files (will be created if not exists).
    
    Optional parameters
    --------------------
    `min_seq_id` : float, optional
        Minimum sequence identity (range 0.0, 1.0)
        - 0.0 (default)

    `s` : float, optional
        Sensitivity.
        - Options: 1.0 (faster), 4.0 (fast), 7.5 (sensitive)
        - Default: 4.0

    `c` : float, optional
        Coverage threshold for alignments
        - 0.8 (default)
        - Determines the minimum fraction of aligned residues required for a match, based on the selected cov_mode

    `cov_mode` : int, optional
        Defines how alignment coverage is calculated:
        - 0: query + target (default)
        - 1: target only
        - 2: query only
        - 3: Target length ≥ x% query length
        - 4: Query length ≥ x% target length
        - 5: Short seq length ≥ x% other seq length
    
    `e` : float, optional
        E-value threshold (range 0.0, inf)
        - 0.001 (default)
    
    `cluster_mode` : int, optional
        Clustering method.
        - 0: Set-Cover (greedy) (default)
        - 1: Connected component (BLASTclust)
        - 2: Greedy by sequence length (CDHIT)

    Output parameters
    ------------------
    `output_mode` : str, optional
        Output mode
        - "to_dict": return a dictionary of clusters
        - "to_list": return a list of clusters
        - "to_gen": return a generator of clusters
        - "to_rel_path": return a list of relative paths to the output files
        - "to_abs_path": return a list of absolute paths to the output files
    
    Returns
    -------
    dict, list, generator, list[str], or None
        - Based on the `output_mode` parameter
    """

    config = EasyClusterConfig(
        fasta_files=fasta_files,
        cluster_prefix=cluster_prefix,
        tmp_dir=tmp_dir,
        min_seq_id=min_seq_id,
        s=s,
        c=c,
        cov_mode=cov_mode,
        e=e,
        cluster_mode=cluster_mode
    )

    config.run()

    return EasyClusterParser(config.cluster_prefix)
