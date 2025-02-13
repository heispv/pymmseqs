# pymmseqs/commands/easy_cluster.py

from pathlib import Path
from typing import Union, Tuple

from pymmseqs import run_mmseqs_command
from pymmseqs.utils import (
    get_caller_dir,
    resolve_path,
    add_arg,
    add_twin_arg
)

def easy_cluster(
    # Required parameters
    input_fasta: Union[str, Path],
    cluster_prefix: Union[str, Path],
    tmp_dir: Union[str, Path],

    # Prefilter parameters
    seed_sub_mat: Tuple[str, str] = ("aa:VTML80.out", "nucl:nucleotide.out"),
    s: float = 4.0,
    k: int = 0,
    target_search_mode: int = 0,
    k_score: Tuple[str, str] = ("seq:2147483647", "prof:2147483647"),
    alph_size: Tuple[str, str] = ("aa:21", "nucl:5"),
    max_seqs: int = 20,
    split: int = 0,
    split_mode: int = 2,
    split_memory_limit: str = "0",
    comp_bias_corr: bool = True,
    comp_bias_corr_scale: float = 1.0,
    diag_score: bool = True,
    exact_kmer_matching: bool = False,
    mask: bool = True,
    mask_prob: float = 0.9,
    mask_lower_case: bool = False,
    min_ungapped_score: int = 15,
    add_self_matches: bool = False,
    spaced_kmer_mode: int = 1,
    spaced_kmer_pattern: str = "",
    local_tmp: Union[str, Path] = "",

    # Alignment parameters
    c: float = 0.8,
    cov_mode: int = 0,
    a: bool = False,
    alignment_mode: int = 3,
    alignment_output_mode: int = 0,
    wrapped_scoring: bool = False,
    e: float = 0.001,
    min_seq_id: float = 0.0,
    min_aln_len: int = 0,
    seq_id_mode: int = 0,
    alt_ali: int = 0,
    max_rejected: int = 2147483647,
    max_accept: int = 2147483647,
    score_bias: float = 0.0,
    realign: bool = False,
    realign_score_bias: float = -0.2,
    realign_max_seqs: int = 2147483647,
    corr_score_weight: float = 0.0,
    gap_open: Tuple[str, str] = ("aa:11", "nucl:5"),
    gap_extend: Tuple[int, int] = ("aa:1", "nucl:2"),
    zdrop: int = 40,

    # Clustering parameters
    cluster_mode: int = 0,
    max_iterations: int = 1000,
    similarity_type: int = 2,
    single_step_clustering: bool = False,
    cluster_steps: int = 3,
    cluster_reassign: bool = False,

    # K-mer matcher parameters
    weights: str = "",
    cluster_weight_threshold: float = 0.9,
    kmer_per_seq: int = 21,
    kmer_per_seq_scale: Tuple[str, str] = ("aa:0.0", "nucl:0.2"),
    adjust_kmer_len: bool = False,
    hash_shift: int = 67,
    include_only_extendable: bool = False,
    ignore_multi_kmer: bool = False,

    # Profile parameters
    pca: float = 0.0,
    pcb: float = 0.0,

    # Misc parameters
    taxon_list: str = "",
    rescore_mode: int = 0,
    dbtype: int = 0,
    shuffle: bool = True,
    createdb_mode: int = 1,

    # Common parameters
    sub_mat: Tuple[str, str] = ("aa:blosum62.out", "nucl:nucleotide.out"),
    max_seq_len: int = 65535,
    db_load_mode: int = 0,
    threads: int = 14,
    compressed: bool = False,
    v: int = 3,
    remove_tmp_files: bool = True,
    force_reuse: bool = False,
    mpi_runner: str = "",


    # Expert parameters
    filter_hits: bool = False,
    sort_results: int = 0,
    write_lookup: bool = False,
) -> None:
    """
    Perform sequence clustering using MMseqs2.

    Parameters
    ----------
    **input_fasta** : Union[str, Path]
        Path to MMseqs2 sequence database created with createdb.

    **cluster_prefix** : Union[str, Path]
        Output cluster database path prefix (will create multiple files with this prefix).

    **tmp_dir** : Union[str, Path]
        Temporary directory for intermediate files (will be created if not exists).

    Prefilter Parameters
    --------------------
    **seed_sub_mat** : Tuple[str, str], optional
        Substitution matrix for k-mer generation as (type:path, type:path)
        - Default: ("aa:VTML80.out", "nucl:nucleotide.out")

    **s** : float, optional
        Sensitivity.
        - Options: 1.0 (faster), 4.0 (fast), 7.5 (sensitive)
        - Default: 4.0

    **k** : int, optional
        k-mer length.
        - 0: automatically set to optimum (default)

    **target_search_mode** : int, optional
        Target search mode.
        - 0: regular k-mer (default)
        - 1: similar k-mer

    **k_score** : Tuple[str, str], optional
        k-mer thresholds for sequence and profile searches.
        - Default: ("seq:2147483647", "prof:2147483647")

    **alph_size** : Tuple[str, str], optional
        Alphabet sizes for amino acid (protein) and nucleotide sequences (range 2-21)
        - ("aa:21", "nucl:5") (default)
            - aa:21: 20 amino acids + X for unknown residues
            - nucl:5: 4 nucleotides + N for unknown bases

    **max_seqs** : int, optional
        Maximum results per query sequence passing prefilter.
        - Default: 20

    **split** : int, optional
        Split input into N chunks
        - 0: set the best split automatically (default)

    **split_mode** : int, optional
        Split strategy
        - 0: split target db
        - 1: split query db
        - 2: auto, depending on main memory (default)

    **split_memory_limit** : str, optional
        Maximum memory allocated per split for processing
        - "0":  all available (default)
            - Use suffixes like K, M, or G (e.g., "4G" for 4 gigabytes)

    **comp_bias_corr** : bool, optional
        Correct for locally biased amino acid composition
        - True (default)
        - False

    **comp_bias_corr_scale** : float, optional
        Scale factor for composition bias correction
        - Range 0, 1
        - 1.0 (default)

    **diag_score** : bool, optional
        Use ungapped diagonal scoring during prefilter
        - True (default)
        - False

    **exact_kmer_matching** : bool, optional
        Extract only exact k-mers for matching
        - True
        - False (default)

    **mask** : bool, optional
        Use low complexity masking
        - True (default)
        - False

    **mask_prob** : float, optional
        Probability threshold for masking low-complexity regions in sequences
        - 0.9 (default)
        - Sequences with low-complexity regions above this threshold are masked during k-mer matching

    **mask_lower_case** : bool, optional
        Mask lowercase letters in k-mer search.
        - True
        - False (default)

    **min_ungapped_score** : int, optional
        Minimum ungapped alignment score
        - 15 (default)
        - Higher values increase specificity but may reduce sensitivity

    **add_self_matches** : bool, optional
        Add entries of queries with themselves for clustering.
        - Default: False

    **spaced_kmer_mode** : int, optional
        Spaced k-mer mode
        - 0: consecutive
        - 1: spaced (default)

    **spaced_kmer_pattern** : str, optional
        Custom pattern for spaced k-mers used during k-mer matching.
        - Define a pattern of 1s (match positions) and 0s (ignore positions)
        - Example: "1101011" means 5 match positions and 2 ignored positions
        - Increases sensitivity by focusing on conserved regions while allowing flexibility in less conserved areas.

    **local_tmp** : str, optional
        Path to an alternative temporary directory for storing intermediate files
        - Useful for reducing I/O load on shared storage systems (e.g., NFS)
        - Default: Temporary files are stored in the main tmpDir


    Alignment Parameters
    --------------------
    c** : float, optional
        Coverage threshold for alignments
        - 0.8 (default)
        - Determines the minimum fraction of aligned residues required for a match, based on the selected cov_mode

    **cov_mode** : int, optional
        Defines how alignment coverage is calculated:
        - 0: query + target (default)
        - 1: target only
        - 2: query only
        - 3: Target length ≥ x% query length
        - 4: Query length ≥ x% target length
        - 5: Short seq length ≥ x% other seq length

    a** : bool, optional
        Add backtrace string (convert to alignments with mmseqs convertalis module)
        - True
        - False (default)

    **alignment_mode** : int, optional
        Alignment detail level
        - 0: auto
        - 1: score + end_po
        - 2: + start_pos + cov
        - 3: + seq.id (default)
        - 4: only ungapped alignment

    **alignment_output_mode** : int, optional
        Output detail level
        - 0: auto (default)
        - 1: score + end_po
        - 2: + start_pos + cov
        - 3: + seq.id
        - 4: only ungapped alignment
        - 5: score only (output) cluster format

    **wrapped_scoring** : bool, optional
        Enable wrapped diagonal scoring for nucleotide sequences by doubling the query sequence
        - True
        - False (default)

    **e** : float, optional
        E-value threshold (range 0.0, inf)
        - 0.001 (default)

    **min_seq_id** : float, optional
        Minimum sequence identity (range 0.0, 1.0)
        - 0.0 (default)

    **min_aln_len** : int, optional
        Minimum alignment length (range 0, inf)
        - 0 (default)

    **seq_id_mode** : int, optional
        Defines how sequence identity calculation is based on
        - 0: Alignment length (default)
        - 1: Shorter sequence
        - 2: Longer sequence

    **alt_ali** : int, optional
        Number of alternative alignments to show
        - 0 (default)

    **max_rejected** : int, optional
        Maximum rejected alignments before alignment calculation for a query is stopped
        - 2147483647 (default)

    **max_accept** : int, optional
        Maximum accepted alignments before alignment calculation for a query is stopped
        - 2147483647 (default)

    **score_bias** : float, optional
        Score bias added to alignment scores (in bits)
        - 0.0: no bias (default)
        - Adjusts alignment scores to favor or penalize certain alignments

    **realign** : bool, optional
        Compute more conservative, shorter alignments (scores and E-values not changed)
        - True
        - False (default)

    **realign_score_bias** : float, optional
        Additional score bias applied during realignment to compute more conservative alignments
        - -0.2 (default)
        - A negative value encourages shorter, more precise alignments

    **realign_max_seqs** : int, optional
        Maximum number of results to return in realignment
        - 2147483647 (default)

    **ccorr_score_weight** : float, optional
        Weight of backtrace correlation score that is added to the alignment score
        - 0.0 (default)
        - Higher values increase the influence of the backtrace correlation on the final alignment score

    **gap_open** : Tuple[str, str], optional
        Gap open costs for amino acid (protein) and nucleotide alignments
        - ("aa:11", "nucl:5") (default)
            - aa:x: Gap open cost for protein alignments
            - nucl:x: Gap open cost for nucleotide alignments
        - Higher values penalize gap openings more heavily, favoring fewer gaps in alignments

    **gap_extend** : Tuple[str, str], optional
        Gap extension costs for amino acid (protein) and nucleotide alignments
        - ("aa:1", "nucl:2") (default)
            - aa:x: Cost for extending a gap in protein alignments
            - nucl:x: Cost for extending a gap in nucleotide alignments
        - Lower values allow longer gaps; higher values penalize gaps more heavily

    **zdrop** : int, optional
        Maximum score drop allowed before truncating the alignment (nucleotide alignments only)
        - 40 (default)
            - Terminates alignments early in low-quality regions to improve computational efficiency

    Clustering Parameters
    ----------------------
    **cluster_mode** : int, optional
        Clustering method.
        - 0: Set-Cover (greedy) (default)
        - 1: Connected component (BLASTclust)
        - 2: Greedy by sequence length (CDHIT)

    **max_iterations** : int, optional
        Maximum depth of breadth first search for connected component clustering.
        - 1000 (default)

    **similarity_type** : int, optional
        Type of score used for clustering.
        - 1: Alignment score
        - 2: Sequence identity (default)

    **single_step_clustering** : bool, optional
        Use simple clustering without cascading.
        - True
        - False: Uses cascading (default)

    **cluster_steps** : int, optional
        Number of clustering steps for cascaded clustering (1 to `s` option).
        - 3 (default)

    **cluster_reassign** : bool, optional
        Cascaded clustering can cluster sequence that do not fulfill the clustering criteria.
        Enable cluster reassignment to correct errors.
        - True
        - False (default)

    K-mer Matcher Parameters
    -------------------------
    **weights** : str, optional
        File containing sequence weights for cluster prioritization.
        Each line in the file corresponds to a sequence in the input
        file, with a floating-point weight determining its influence on clustering (higher values = higher priority)
        - Default: ""
        - Example: A text file like this:\n1.0\n0.5\n0.2\n

    **cluster_weight_threshold** : float, optional
        Weight threshold for cluster prioritization.
        - 0.9 (default)

    **kmer_per_seq** : int, optional
        Number of k-mers per sequence.
        - 21 (default)

    **kmer_per_seq_scale** : Tuple[str, str], optional
        Scale k-mer per sequence based on sequence length as (kmer-pers-seq val + scale * seq-len)
        - ("aa:0.0", "nucl:0.2") (default)

    **adjust_kmer_len** : bool, optional
        Adjust k-mer length based on specificity (only for nucleotide).
        - True
        - False (default)

    **hash_shift** : int, optional
        Shift k-mer hash initialization.
        - 67 (default)

    **include_only_extendable** : bool, optional
        Include only extendable k-mers.
        - True
        - False (default)

    **ignore_multi_kmer** : bool, optional
        Skip k-mers occurring multiple times (>=2)
        - True
        - False (default)

    Profile Parameters
    -----------------
    **pca** : float, optional
        Pseudo count admixture strength for profile construction
        - 0.0 (default)
        - Higher values increase the weight of pseudo counts, making the profile more conservative
        - Lower values reduce their influence, making the profile more specific to the input sequences

    **pcb** : float, optional
        Controls the threshold for pseudo-count admixture based on the effective number of sequences (Neff) (range 0.0, inf)
        - 0.0 (default)
        - Lower values apply pseudo-counts more aggressively

    Misc Parameters
    ---------------
    **taxon_list** : str, optional
        Taxonomy IDs to filter results by. Multiple IDs can be provided, separated by commas (no spaces)
        - "" (default)
        - Example: "9606,10090"

    **rescore_mode** : int, optional
        Rescore diagonals with:
        - 0: Hamming distance (default)
        - 1: local alignment (score only)
        - 2: local alignment
        - 3: global alignment
        - 4: longest alignment fulfilling window quality criterion
    
    **dbtype** : int, optional  
        Database type
        - 0: Auto-detect (default)
        - 1: Amino acid
        - 2: Nucleotide

    **shuffle** : bool, optional
        Shuffle the input database before processing.
        - True (default)
        - False

    **createdb_mode** : int, optional  
        Database creation mode
        - 0: Copy data
        - 1: Soft link data and write a new index (only works with single-line FASTA/Q files) (default)

    Common Parameters
    -----------------
    **sub_mat** : Tuple[str, str], optional
        Substitution matrix (type:path, type:path)
        type: "aa" or "nucl"
        path: matrix file path
        - ("aa:blosum62.out", "nucl:nucleotide.out")

        Note: find available matrices in the MMseqs2 data directory: (https://github.com/soedinglab/MMseqs2/tree/master/data)

    **max_seq_len** : int, optional
        Maximum sequence length
        - 65535 (default)

    **db_load_mode** : int, optional
        Database preloading method
        - 0: auto (default)
        - 1: fread
        - 2: mmap
        - 3: mmap+touch

    **threads** : int, optional
        CPU threads
        - 14 (default)

    **compressed** : bool, optional
        Compress output
        - True
        - False (default)

    v** : int, optional
        Output verbosity
        - 0: quiet
        - 1: +errors
        - 2: +warnings
        - 3: +info (default)

    **remove_tmp_files** : bool, optional
        Delete temporary files
        - True (default)
        - False

    **force_reuse** : bool, optional
        Reuse tmp filse in tmp/latest folder ignoring parameters and version changes
        - True
        - False (default)
    
    **mpi_runner** : str, optional
        Use MPI on compute cluster with this MPI command (e.g., "mpirun -np 42")
        - "" (default)

    Expert Parameters
    -----------------
    **filter_hits** : bool, optional
        Filter hits by sequence ID and coverage
        - True
        - False (default)

    **sort_results** : int, optional
        Result sorting method
        - 0: No sorting (default)
        - 1: E-value (Alignment) or sequence ID (Hamming) 
    
   **write_lookup** : bool, optional
        Create a `.lookup` file mapping internal IDs to FASTA IDs
        - True (default)
        - False

    Returns
    -------
    None
        Creates the cluster database files at the specified path.

    Raises
    ------
    FileNotFoundError
        If input databases are missing.

    ValueError
        For invalid parameter combinations or values.

    Examples
    --------
    Basic clustering:
    >>> easy_cluster(
        "examples/DB.fasta",
        "results",
        "tmp"
    )
    """

    # TODO: You should add all the checks for the parameters here
    # Validate numerical parameters
    if not (0 <= comp_bias_corr_scale <= 1):
        raise ValueError("comp_bias_corr_scale must be between 0 and 1")
    if not (1.0 <= s <= 7.5):
        raise ValueError("Sensitivity (-s) must be between 1.0 and 7.5")
    if min_seq_id < 0.0 or min_seq_id > 1.0:
        raise ValueError("min_seq_id must be between 0.0 and 1.0")

    # Get the directory of the calling script
    caller_dir = get_caller_dir()

    input_fasta_path = resolve_path(input_fasta, caller_dir)
    cluster_prefix_path = resolve_path(cluster_prefix, caller_dir)
    tmp_dir_path = resolve_path(tmp_dir, caller_dir)

    local_tmp_path = resolve_path(local_tmp, caller_dir) if local_tmp else ""

    # Validate input sequence database
    if not input_fasta_path.exists():
        raise FileNotFoundError(f"Sequence database not found: {input_fasta_path}")

    # Build base command
    args = [
        "easy-cluster",
        str(input_fasta_path),
        str(cluster_prefix_path),
        str(tmp_dir_path),
    ]

    # Add parameters
    # Prefilter
    add_twin_arg(args, "--seed-sub-mat", seed_sub_mat, ("aa:VTML80.out", "nucl:nucleotide.out"), ",")
    add_arg(args, "-s", s, 4.0)
    add_arg(args, "-k", k, 0)
    add_arg(args, "--target-search-mode", target_search_mode, 0)
    add_twin_arg(args, "--k-score", k_score, ("seq:2147483647", "prof:2147483647"), ",")
    add_twin_arg(args, "--alph-size", alph_size, ("aa:21", "nucl:5"), ",")
    add_arg(args, "--max-seqs", max_seqs, 20)
    add_arg(args, "--split", split, 0)
    add_arg(args, "--split-mode", split_mode, 2)
    add_arg(args, "--split-memory-limit", split_memory_limit, "0")
    add_arg(args, "--comp-bias-corr", comp_bias_corr, True)
    add_arg(args, "--comp-bias-corr-scale", comp_bias_corr_scale, 1.0)
    add_arg(args, "--diag-score", diag_score, True)
    add_arg(args, "--exact-kmer-matching", exact_kmer_matching, False)
    add_arg(args, "--mask", mask, True)
    add_arg(args, "--mask-prob", mask_prob, 0.9)
    add_arg(args, "--mask-lower-case", mask_lower_case, False)
    add_arg(args, "--min-ungapped-score", min_ungapped_score, 15)
    add_arg(args, "--add-self-matches", add_self_matches, False)
    add_arg(args, "--spaced-kmer-mode", spaced_kmer_mode, 1)
    add_arg(args, "--spaced-kmer-pattern", spaced_kmer_pattern, "")
    add_arg(args, "--local-tmp", local_tmp_path, "")

    # Alignment
    add_arg(args, "-c", c, 0.8)
    add_arg(args, "--cov-mode", cov_mode, 0)
    add_arg(args, "-a", a, False)
    add_arg(args, "--alignment-mode", alignment_mode, 3)
    add_arg(args, "--alignment-output-mode", alignment_output_mode, 0)
    add_arg(args, "--wrapped-scoring", wrapped_scoring, False)
    add_arg(args, "-e", e, 0.001)
    add_arg(args, "--min-seq-id", min_seq_id, 0.0)
    add_arg(args, "--min-aln-len", min_aln_len, 0)
    add_arg(args, "--seq-id-mode", seq_id_mode, 0)
    add_arg(args, "--alt-ali", alt_ali, 0)
    add_arg(args, "--max-rejected", max_rejected, 2147483647)
    add_arg(args, "--max-accept", max_accept, 2147483647)
    add_arg(args, "--score-bias", score_bias, 0.0)
    add_arg(args, "--realign", realign, False)
    add_arg(args, "--realign-score-bias", realign_score_bias, -0.2)
    add_arg(args, "--realign-max-seqs", realign_max_seqs, 2147483647)
    add_arg(args, "--corr-score-weight", corr_score_weight, 0.0)
    add_twin_arg(args, "--gap-open", gap_open, ("aa:11", "nucl:5"), ",")
    add_twin_arg(args, "--gap-extend", gap_extend, ("aa:1", "nucl:2"), ",")
    add_arg(args, "--zdrop", zdrop, 40)

    # Clustering
    add_arg(args, "--cluster-mode", cluster_mode, 0)
    add_arg(args, "--max-iterations", max_iterations, 1000)
    add_arg(args, "--similarity-type", similarity_type, 2)
    add_arg(args, "--single-step-clustering", single_step_clustering, False)
    add_arg(args, "--cluster-steps", cluster_steps, 3)
    add_arg(args, "--cluster-reassign", cluster_reassign, False)

    # K-mer Matcher
    add_arg(args, "--weights", weights, "")
    add_arg(args, "--cluster-weight-threshold", cluster_weight_threshold, 0.9)
    add_arg(args, "--kmer-per-seq", kmer_per_seq, 21)
    add_twin_arg(args, "--kmer-per-seq-scale", kmer_per_seq_scale, ("aa:0.0", "nucl:0.2"), ",")
    add_arg(args, "--adjust-kmer-len", adjust_kmer_len, False)
    add_arg(args, "--hash-shift", hash_shift, 67)
    add_arg(args, "--include-only-extendable", include_only_extendable, False)
    add_arg(args, "--ignore-multi-kmer", ignore_multi_kmer, False)

    # Profile
    add_arg(args, "--pca", pca, 0.0)
    add_arg(args, "--pcb", pcb, 0.0)

    # Misc
    add_arg(args, "--taxon-list", taxon_list, "")
    add_arg(args, "--rescore-mode", rescore_mode, 0)
    add_arg(args, "--dbtype", dbtype, 0)
    add_arg(args, "--shuffle", shuffle, True)
    add_arg(args, "--createdb-mode", createdb_mode, 1)

    # Common
    add_twin_arg(args, "--sub-mat", sub_mat, ("aa:blosum62.out", "nucl:nucleotide.out"), ",")
    add_arg(args, "--max-seq-len", max_seq_len, 65535)
    add_arg(args, "--db-load-mode", db_load_mode, 0)
    add_arg(args, "--threads", threads, 14)
    add_arg(args, "--compressed", compressed, False)
    add_arg(args, "-v", v, 3)
    add_arg(args, "--remove-tmp-files", remove_tmp_files, True)
    add_arg(args, "--force-reuse", force_reuse, False)
    add_arg(args, "--mpi-runner", mpi_runner, "")

    # Expert
    add_arg(args, "--filter-hits", filter_hits, False)
    add_arg(args, "--sort-results", sort_results, 0)
    add_arg(args, "--write-lookup", write_lookup, False)

    # Execute command
    mmseqs_output = run_mmseqs_command(args)
    print(mmseqs_output.stdout)
    if mmseqs_output.stderr:
        print(mmseqs_output.stderr)
    print(f"MMseqs2 clustering completed. Results saved to: {cluster_prefix_path}")
