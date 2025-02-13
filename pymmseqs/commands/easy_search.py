# pymmseqs/commands/easy_search.py

from pathlib import Path
from typing import Union, List, Tuple

from pymmseqs import run_mmseqs_command
from pymmseqs.utils import (
    get_caller_dir,
    resolve_path,
    add_arg,
    add_twin_arg
)


def easy_search(
    # Required parameters
    query_fasta: Union[str, Path],
    target_fasta_or_db: Union[str, Path],
    alignment_file: Union[str, Path],
    tmp_dir: Union[str, Path],
    
    # Prefilter parameters
    comp_bias_corr: bool = True,
    comp_bias_corr_scale: float = 1.0,
    add_self_matches: bool = False,
    seed_sub_mat: Tuple[str, str] = ("aa:VTML80.out", "nucl:nucleotide.out"),
    s: float = 5.7,
    k: int = 0,
    target_search_mode: int = 0,
    k_score: Tuple[str, str] = ("seq:2147483647", "prof:2147483647"),
    alph_size: Tuple[str, str] = ("aa:21", "nucl:5"),
    max_seqs: int = 300,
    split: int = 0,
    split_mode: int = 2,
    split_memory_limit: str = "0",
    diag_score: bool = True,
    exact_kmer_matching: bool = False,
    mask: bool = True,
    mask_prob: float = 0.9,
    mask_lower_case: bool = False,
    min_ungapped_score: int = 15,
    spaced_kmer_mode: int = 1,
    spaced_kmer_pattern: str = "",
    local_tmp: Union[str, Path] = "",
    disk_space_limit: str = "0",
    
    # Alignment parameters
    a: bool = False,
    alignment_mode: int = 3,
    alignment_output_mode: int = 0,
    wrapped_scoring: bool = False,
    e: float = 0.001,
    min_seq_id: float = 0.0,
    min_aln_len: int = 0,
    seq_id_mode: int = 0,
    alt_ali: int = 0,
    c: float = 0.0,
    cov_mode: int = 0,
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
    exhaustive_search_filter: bool = False,
    
    # Profile parameters
    pca: float = 0.0,
    pcb: float = 0.0,
    mask_profile: bool = True,
    e_profile: float = 0.1,
    wg: bool = False,
    filter_msa: bool = True,
    filter_min_enable: int = 0,
    profile_max_seq_id: float = 0.9,
    qid: str = "0.0",
    qsc: float = -20.0,
    profile_cov: float = 0.0,
    diff: int = 1000,
    pseudo_cnt_mode: int = 0,
    num_iterations: int = 1,
    exhaustive_search: bool = False,
    lca_search: bool = False,
    
    # Misc parameters
    taxon_list: str = "",
    prefilter_mode: int = 0,
    rescore_mode: int = 0,
    allow_deletion: bool = False,
    min_length: int = 30,
    max_length: int = 32734,
    max_gaps: int = 2147483647,
    contig_start_mode: int = 2,
    contig_end_mode: int = 2,
    orf_start_mode: int = 1,
    forward_frames: List[int] = [1, 2, 3],
    reverse_frames: List[int] = [1, 2, 3],
    translation_table: int = 1,
    translate: bool = False,
    use_all_table_starts: bool = False,
    id_offset: int = 0,
    sequence_overlap: int = 0,
    sequence_split_mode: int = 1,
    headers_split_mode: int = 0,
    search_type: int = 0,
    start_sens: float = 4.0,
    sens_steps: int = 1,
    translation_mode: int = 0,
    format_mode: int = 0,
    format_output: List[str] = [
        "query",
        "target",
        "fident",
        "alnlen",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "evalue",
        "bits"
    ],
    overlap: float = 0.0,
    dbtype: int = 0,
    shuffle: bool = True,
    createdb_mode: int = 0,
    greedy_best_hits: bool = False,
    
    # Common parameters
    sub_mat: Tuple[str, str] = ("aa:blosum62.out", "nucl:nucleotide.out"),
    max_seq_len: int = 65535,
    db_load_mode: int = 0,
    threads: int = 14,
    compressed: bool = False,
    v: int = 3,
    gpu: bool = False,
    gpu_server: bool = False,
    mpi_runner: str = "",
    force_reuse: bool = False,
    remove_tmp_files: bool = True,
    
    # Expert parameters
    filter_hits: bool = False,
    sort_results: int = 0,
    create_lookup: bool = False,
    chain_alignments: bool = False,
    merge_query: bool = True,
    strand: int = 1,
    db_output: bool = False,
    write_lookup: bool = False,
) -> None:
    """
    Perform a fast and sensitive sequence search using MMseqs2 easy-search.

    Parameters
    ----------
    **query_fasta** : Union[str, Path]
        Path to one or more query FASTA files. Can be compressed with .gz or .bz2.

    **target_fasta_or_db** : Union[str, Path]
        Path to a target FASTA file (optionally compressed) or an MMseqs2 target database.

    **alignment_output** : Union[str, Path]
        Path to the output file where alignments will be stored.

    **tmp_dir** : Union[str, Path]
        Temporary directory for intermediate files. Will be created if not existing.

    Prefilter Parameters
    -------------------
    **comp_bias_corr** : bool, optional
        Correct for locally biased amino acid composition
        - True (default)
        - False

    **comp_bias_corr_scale** : float, optional
        Scale factor for composition bias correction
        - Range 0, 1
        - 1.0 (default)

    **add_self_matches** : bool, optional
        Artificially add entries of queries with themselves (for clustering)
        - True
        - False (default)

    **seed_sub_mat** : Tuple[str, str], optional
        Substitution matrix for k-mer generation as (type:path, type:path)
        type: "aa" or "nucl"
        path: matrix file path
        - ("aa:VTML80.out", "nucl:nucleotide.out") (default)

        Note: find available matrices in the MMseqs2 data directory: (https://github.com/soedinglab/MMseqs2/tree/master/data)

    **s** : float, optional
        Sensitivity
        - 1.0: faster
        - 4.0: fast
        - 5.7 (default)
        - 7.5: sensitive

    **k** : int, optional
        k-mer length
        - 0: auto (default)

    **target_search_mode** : int, optional
        Target search mode
        - 0: regular k-mer (default)
        - 1: similar k-mer

    **k_score** : Tuple[str, str], optional
        k-mer thresholds for sequence and profile searches
        - ("seq:2147483647", "prof:2147483647") (default)

    **alph_size** : Tuple[str, str], optional
        Alphabet sizes for amino acid (protein) and nucleotide sequences (range 2-21)
        - ("aa:21", "nucl:5") (default)
            - aa:21: 20 amino acids + X for unknown residues
            - nucl:5: 4 nucleotides + N for unknown bases

    **max_seqs** : int, optional
        Maximum results per query passing prefilter
        - 300 (default)
        - Higher values increase sensitivity but may slow down the search

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

    **disk_space_limit** : str, optional
        Max disk space usage
        - "0": unlimited (default)
        - Use suffixes like K, M, or G (e.g., "100G" for 100 gigabytes)

    Alignment Parameters
    --------------------
    **a** : bool, optional
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

    **c** : float, optional
        Coverage threshold for alignments
        - 0.0 (default)
        - Determines the minimum fraction of aligned residues required for a match, based on the selected cov_mode
    
    **cov_mode** : int, optional
        Defines how alignment coverage is calculated:
        - 0: query + target (default)
        - 1: target only
        - 2: query only
        - 3: Target length ≥ x% query length
        - 4: Query length ≥ x% target length
        - 5: Short seq length ≥ x% other seq length

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

    **corr_score_weight** : float, optional
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

    **exhaustive_search_filter** : bool, optional
        Filter result during search
        - True
        - False (default)

    Profile Parameters
    ------------------
    **pca** : float, optional
        Pseudo count admixture strength for profile construction
        - 0.0 (default)
        - Higher values increase the weight of pseudo counts, making the profile more conservative
        - Lower values reduce their influence, making the profile more specific to the input sequences

    **pcb** : float, optional
        Controls the threshold for pseudo-count admixture based on the effective number of sequences (Neff) (range 0.0, inf)
        - 0.0 (default)
        - Lower values apply pseudo-counts more aggressively

    **mask_profile** : bool, optional
        Mask low-complexity regions in the query sequence of a profile using TANTAN
        - True (default)
        - False

    **e_profile** : float, optional
        E-value threshold for including sequence matches in the profile
        - 0.1 (default)

    **wg** : bool, optional
        Use global sequence weighting for profile calculation
        - True
        - False (default)

    **filter_msa** : bool, optional
        Filter MSA
        - True (default)
        - False

    **filter_min_enable** : int, optional
        Minimum number of sequences required to trigger filtering of MSAs
        - 0: Always filter (default)
        - N > 0: Filter only if the MSA contains more than N sequences

    **max_seq_id** : float, optional
        Maximum pairwise sequence identity for redundancy reduction in the output MSA (range 0.0, 1.0)
        - 0.9 (default)
        - Filters sequences to ensure no two sequences in the output share more than the specified identity

    **qid** : str, optional
        Filters output MSAs by minimum sequence identity with the query (range 0.0, 1.0)
        - 0.0: no filtering (default)
        - Can specify multiple thresholds as a comma-separated list (e.g., "0.15,0.30,0.50") to create filter buckets
            - Example: "0.15,0.30,0.50" creates buckets for sequences with identities in ]0.15-0.30] and ]0.30-0.50]

    **qsc** : float, optional
        Filters output MSAs by minimum score per aligned residue with query sequences (range -50.0, 100.0)
        - -20.0 (default)
        - Higher values reduce diversity in the output MSAs by retaining only high-scoring alignments

    **profile_cov** : float, optional
        Minimum fraction of query residues covered by matched sequences to filter output MSAs (range 0.0, 1.0)
        - 0.0 (default)

    **diff** : int, optional
        Filters MSAs by selecting the most diverse sequences, ensuring at least this many sequences are kept in each MSA block of length 50
        - 1000 (default)

    **pseudo_cnt_mode** : int, optional
        Pseudocount method
        - 0: substitution-matrix (default)
        - 1: context-specific pseudocounts

    **num_iterations** : int, optional
        Number of iterative profile search iterations
        - 1: (default)

    **exhaustive_search** : bool, optional
        For bigger profile DB, run iteratively the search by greedily swapping the search results
        - True
        - False (default)

    **lca_search** : bool, optional
        Enable LCA candidate search
        - True
        - False (default)

    Misc Parameters
    ---------------
    **taxon_list** : str, optional
        Taxonomy IDs to filter results by. Multiple IDs can be provided, separated by commas (no spaces)
        - "" (default)
        - Example: "9606,10090"

    **prefilter_mode** : int, optional
        Prefilter method
        - 0: kmer/ungapped (default)
        - 1: ungapped
        - 2: nofilter
        - 3: ungapped+gapped

    **rescore_mode** : int, optional
        Rescore diagonals with:
        - 0: Hamming distance (default)
        - 1: local alignment (score only)
        - 2: local alignment
        - 3: global alignment
        - 4: longest alignment fulfilling window quality criterion

    **allow_deletion** : bool, optional
        Allow deletions in MSA
        - True
        - False (default)

    **min_length** : int, optional
        Minimum codon number in open reading frames (ORFs)
        - 30 (default)

    **max_length** : int, optional
        Maximum codon number in open reading frames (ORFs)
        - 32734 (default)

    **max_gaps** : int, optional
        Maximum number of codons with gaps or unknown residues before an open reading frame is rejected
        - 2147483647 (default)

    **contig_start_mode** : int, optional
        Contig start handling
        - 0: incomplete
        - 1: complete
        - 2: both (default)

    **contig_end_mode** : int, optional
        Contig end handling
        - 0: incomplete
        - 1: complete
        - 2: both (default)

    **orf_start_mode** : int, optional
        ORF start handling
        - 0: from start to stop
        - 1: from any to stop (default)
        - 2: from last encountered start to stop (no start in the middle)

    **forward_frames** : List[int], optional
        Comma-separated list of frames on the forward strand to be extracted
        - [1, 2, 3] (default)

    **reverse_frames** : List[int], optional
        Comma-separated list of frames on the reverse strand to be extracted
        - [1, 2, 3] (default)

    **translation_table** : int, optional  
        Specifies the genetic code table to use 
        - 1: Canonical (default)
        - 2: Vert Mitochondrial
        - 3: Yeast Mitochondrial
        - 4: Mold Mitochondrial
        - 5: Invert Mitochondrial
        - 6: Ciliate
        - 9: Flatworm Mitochondrial
        - 10: Euplotid
        - 11: Prokaryote
        - 12: Alt Yeast
        - 13: Ascidian Mitochondrial
        - 14: Alt Flatworm Mitochondrial
        - 15: Blepharisma
        - 16: Chlorophycean Mitochondrial
        - 21: Trematode Mitochondrial
        - 22: Scenedesmus Mitochondrial
        - 23: Thraustochytrium Mitochondrial
        - 24: Pterobranchia Mitochondrial
        - 25: Gracilibacteria
        - 26: Pachysolen
        - 27: Karyorelict
        - 28: Condylostoma
        - 29: Mesodinium
        - 30: Pertrich
        - 31: Blastocrithidia

    **translate** : bool, optional
        Translate open reading frames (ORFs) to amino acid
        - True
        - False (default)

    **use_all_table_starts** : bool, optional
        Use all start codons
        - True
        - False: only ATG (AUG) (default)

    **id_offset** : int, optional
        Numeric IDs in index file are offset by this value
        - 0 (default)

    **sequence_overlap** : int, optional
        Overlap between sequences
        - 0 (default)

    **sequence_split_mode** : int, optional
        Method for splitting sequences during processing
        - 0: Copy data (creates a full copy of the sequence data).
        - 1: Soft link data and write a new index (saves disk space by linking to the original data) (default)

    **headers_split_mode** : int, optional
        Header split method
        - 0: Split positions (Headers are split based on predefined positions) (default)
        - 1: Original header (Headers are preserved as-is without splitting)

    **search_type** : int, optional
        Search mode:
        - 0: auto (default)
        - 1: amino
        - 2: translated
        - 3: nucleotide
        - 4: translated alignment

    **start_sens** : float, optional
        Initial sensitivity
        - 4.0 (default)

    **sens_steps** : int, optional
        Number of search steps performed from `start_sens` argument to `s` argument
        - 1 (default)

    **translation_mode** : int, optional
        Translation AA seq from nucletoide method
        - 0: Open Reading Frames (ORFs) (default)
        - 1: Full Reading Frames
    
    **format_mode** : int, optional  
        Output format type  
        - 0: BLAST-TAB (default)  
        - 1: SAM  
        - 2: BLAST-TAB + query/db length  
        - 3: Pretty HTML  
        - 4: BLAST-TAB + column headers  

        Notes:  
        - BLAST-TAB (0) and BLAST-TAB + column headers (4) support custom output formats via `format_output`.  

    **format_output** : str, optional  
        Comma-separated list of output columns to include in results.  
        Available columns:  
        - query, target, evalue, gapopen, pident, fident, nident, qstart, qend, qlen  
        - tstart, tend, tlen, alnlen, raw, bits, cigar, qseq, tseq, qheader, theader, qaln, taln  
        - qframe, tframe, mismatch, qcov, tcov, qset, qsetid, tset, tsetid, taxid, taxname, taxlineage  
        - qorfstart, qorfend, torfstart, torfend, ppos  

        - Default: ["query", "target", "fident", "alnlen", "mismatch", "gapopen",
                    "qstart", "qend", "tstart", "tend", "evalue", "bits"]

    **overlap** : float, optional  
        Maximum allowed overlap of covered regions between query and target sequences.  
        - 0.0 (default)  
        - Higher values allow more overlap, reducing redundancy.  

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
        - 0: Copy data (default)  
        - 1: Soft link data and write a new index (only works with single-line FASTA/Q files)  

        Notes:  
        - Using soft links reduces disk usage by referencing the original sequence data instead of duplicating it.  

    **greedy_best_hits** : bool, optional  
        Select best hits greedily to maximize query coverage.  
        - True  
        - False (default)
    
    Common Parameters
    ----------------
    **sub_mat** : Tuple[str, str], optional
        Substitution matrix (type:path, type:path)
        type: "aa" or "nucl"
        path: matrix file path
        - Default: ("aa:blosum62.out", "nucl:nucleotide.out")

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

    **v** : int, optional
        Output verbosity
        - 0: quiet
        - 1: +errors
        - 2: +warnings
        - 3: +info (default)

    **gpu** : bool, optional
        Use GPU (CUDA) if possible
        - True
        - False (default)

    **gpu_server** : bool, optional
        Use GPU server
        - True
        - False (default)

    **mpi_runner** : str, optional
        Use MPI on compute cluster with this MPI command (e.g., "mpirun -np 42")
        - "" (default)

    **force_reuse** : bool, optional
        Reuse tmp filse in tmp/latest folder ignoring parameters and version changes
        - True
        - False (default)

    **remove_tmp_files** : bool, optional
        Delete temporary files
        - True (default)
        - False

    Expert Parameters
    ----------------
    **filter_hits** : bool, optional
        Filter hits by sequence ID and coverage
        - True
        - False (default)

    **sort_results** : int, optional
        Result sorting method
        - 0: No sorting (default)
        - 1: E-value (Alignment) or sequence ID (Hamming) 

    **create_lookup** : bool, optional
        Create lookup file (can be very large)
        - True
        - False (default)

    **chain_alignments** : bool, optional
        Chain overlapping alignments
        - True
        - False (default)

    **merge_query** : bool, optional
        Combine ORFs/split sequences to a single entry
        - True (default)
        - False

    **strand** : int, optional
        Strand selection (only works for DNA/DNA search)
        - 0: reverse
        - 1: forward (default)
        - 2: both
    
    **db_output** : bool, optional
        Return the result as DB instead of a text file
        - True: DB output
        - False: Text file output (default)
    
    **write_lookup** : bool, optional
        Create a `.lookup` file mapping internal IDs to FASTA IDs
        - True (default)
        - False

    Returns
    -------
    None
        Writes alignment output to specified file.

    Raises
    ------
    FileNotFoundError
        If input FASTA files or databases are missing.

    ValueError
        For invalid parameter combinations.

    Examples
    --------
    Basic protein sequence search:
    >>> easy_search(
        query_fasta="query.fasta",
        target_fasta_or_db="target.fasta",
        alignment_output="output.m8",
        tmp_dir="tmp_search",
        s=5.7,
        threads=8
    )

    Iterative profile search:
    >>> easy_search(
        query_fasta="query.fasta",
        target_fasta_or_db="database.fasta",
        alignment_output="profile_search.m8",
        tmp_dir="tmp_search",
        num_iterations=3
    )

    Search against a small profile database:
    >>> easy_search(
        query_fasta="query.fasta",
        target_fasta_or_db="pfam_db",
        alignment_output="pfam_results.m8",
        tmp_dir="tmp_search"
    )

    """
    # TODO: This should be in a file like config check or something
    # Validate numerical parameters
    # if not (0 <= comp_bias_corr_scale <= 1):
    #     raise ValueError("comp_bias_corr_scale must be between 0 and 1")
    # if not (1.0 <= s <= 7.5):
    #     raise ValueError("Sensitivity (-s) must be between 1.0 and 7.5")
    # if min_seq_id < 0.0 or min_seq_id > 1.0:
    #     raise ValueError("min_seq_id must be between 0.0 and 1.0")

    # Get the directory of the calling script
    caller_dir = get_caller_dir()

    query_fasta_path = resolve_path(query_fasta, caller_dir)
    target_fasta_or_db_path = resolve_path(target_fasta_or_db, caller_dir)
    alignment_file_path = resolve_path(alignment_file, caller_dir)
    tmp_dir_path = resolve_path(tmp_dir, caller_dir)

    local_tmp_path = resolve_path(local_tmp, caller_dir) if local_tmp else ""

    # Validate input files
    for fasta_path, name in [(query_fasta_path, "Query"), (target_fasta_or_db_path, "Target")]:
        if not fasta_path.exists():
            raise FileNotFoundError(f"{name} FASTA file not found: {fasta_path}")

    # Build base command
    args = [
        "easy-search",
        str(query_fasta_path),
        str(target_fasta_or_db_path),
        str(alignment_file_path),
        str(tmp_dir_path),
    ]

    # Prefilter
    add_arg(args, "--comp-bias-corr", comp_bias_corr, True)
    add_arg(args, "--comp-bias-corr-scale", comp_bias_corr_scale, 1.0)
    add_arg(args, "--add-self-matches", add_self_matches, False)
    add_twin_arg(args, "--seed-sub-mat", seed_sub_mat, ("aa:VTML80.out", "nucl:nucleotide.out"), ",")
    add_arg(args, "-s", s, 5.7)
    add_arg(args, "-k", k, 0)
    add_arg(args, "--target-search-mode", target_search_mode, 0)
    add_twin_arg(args, "--k-score", k_score, ("seq:2147483647", "prof:2147483647"), ",")
    add_twin_arg(args, "--alph-size", alph_size, ("aa:21", "nucl:5"), ",")
    add_arg(args, "--max-seqs", max_seqs, 300)
    add_arg(args, "--split", split, 0)
    add_arg(args, "--split-mode", split_mode, 2)
    add_arg(args, "--split-memory-limit", split_memory_limit, "0")
    add_arg(args, "--diag-score", diag_score, True)
    add_arg(args, "--exact-kmer-matching", exact_kmer_matching, False)
    add_arg(args, "--mask", mask, True)
    add_arg(args, "--mask-prob", mask_prob, 0.9)
    add_arg(args, "--mask-lower-case", mask_lower_case, False)
    add_arg(args, "--min-ungapped-score", min_ungapped_score, 15)
    add_arg(args, "--spaced-kmer-mode", spaced_kmer_mode, 1)
    add_arg(args, "--spaced-kmer-pattern", spaced_kmer_pattern, "")
    add_arg(args, "--local-tmp", local_tmp_path, "")
    add_arg(args, "--disk-space-limit", disk_space_limit, "0")

    # Alignment
    add_arg(args, "-a", a, False)
    add_arg(args, "--alignment-mode", alignment_mode, 3)
    add_arg(args, "--alignment-output-mode", alignment_output_mode, 0)
    add_arg(args, "--wrapped-scoring", wrapped_scoring, False)
    add_arg(args, "-e", e, 0.001)
    add_arg(args, "--min-seq-id", min_seq_id, 0.0)
    add_arg(args, "--min-aln-len", min_aln_len, 0)
    add_arg(args, "--seq-id-mode", seq_id_mode, 0)
    add_arg(args, "--alt-ali", alt_ali, 0)
    add_arg(args, "-c", c, 0.0)
    add_arg(args, "--cov-mode", cov_mode, 0)
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
    add_arg(args, "--exhaustive-search-filter", exhaustive_search_filter, False)

    # Profile
    add_arg(args, "--pca", pca, 0.0)
    add_arg(args, "--pcb", pcb, 0.0)
    add_arg(args, "--mask-profile", mask_profile, True)
    add_arg(args, "--e-profile", e_profile, 0.1)
    add_arg(args, "--wg", wg, False)
    add_arg(args, "--filter-msa", filter_msa, True)
    add_arg(args, "--filter-min-enable", filter_min_enable, 0)
    add_arg(args, "--max-seq-id", profile_max_seq_id, 0.9)
    add_arg(args, "--qid", qid, "0.0")
    add_arg(args, "--qsc", qsc, -20.0)
    add_arg(args, "--cov", profile_cov, 0.0)
    add_arg(args, "--diff", diff, 1000)
    add_arg(args, "--pseudo-cnt-mode", pseudo_cnt_mode, 0)
    add_arg(args, "--num-iterations", num_iterations, 1)
    add_arg(args, "--exhaustive-search", exhaustive_search, False)
    add_arg(args, "--lca-search", lca_search, False)

    # Misc
    add_arg(args, "--taxon-list", taxon_list, "")
    add_arg(args, "--prefilter-mode", prefilter_mode, 0)
    add_arg(args, "--rescore-mode", rescore_mode, 0)
    add_arg(args, "--allow-deletion", allow_deletion, False)
    add_arg(args, "--min-length", min_length, 30)
    add_arg(args, "--max-length", max_length, 32734)
    add_arg(args, "--max-gaps", max_gaps, 2147483647)
    add_arg(args, "--contig-start-mode", contig_start_mode, 2)
    add_arg(args, "--contig-end-mode", contig_end_mode, 2)
    add_arg(args, "--orf-start-mode", orf_start_mode, 1)
    add_arg(args, "--forward-frames", ",".join(map(str, forward_frames)), "1,2,3")
    add_arg(args, "--reverse-frames", ",".join(map(str, reverse_frames)), "1,2,3")
    add_arg(args, "--translation-table", translation_table, 1)
    add_arg(args, "--translate", translate, False)
    add_arg(args, "--use-all-table-starts", use_all_table_starts, False)
    add_arg(args, "--id-offset", id_offset, 0)
    add_arg(args, "--sequence-overlap", sequence_overlap, 0)
    add_arg(args, "--sequence-split-mode", sequence_split_mode, 1)
    add_arg(args, "--headers-split-mode", headers_split_mode, 0)
    add_arg(args, "--search-type", search_type, 0)
    add_arg(args, "--start-sens", start_sens, 4.0)
    add_arg(args, "--sens-steps", sens_steps, 1)
    add_arg(args, "--translation-mode", translation_mode, 0)
    add_arg(args, "--format-mode", format_mode, 0)
    add_arg(
    args, "--format-output", ",".join(map(str, format_output)),
    "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
    )
    add_arg(args, "--overlap", overlap, 0.0)
    add_arg(args, "--dbtype", dbtype, 0)
    add_arg(args, "--shuffle", shuffle, True)
    add_arg(args, "--createdb-mode", createdb_mode, 0)
    add_arg(args, "--greedy-best-hits", greedy_best_hits, False)

    # Common
    add_twin_arg(args, "--sub-mat", sub_mat, ("aa:blosum62.out", "nucl:nucleotide.out"), ",")
    add_arg(args, "--max-seq-len", max_seq_len, 65535)
    add_arg(args, "--db-load-mode", db_load_mode, 0)
    add_arg(args, "--threads", threads, 14)
    add_arg(args, "--compressed", compressed, False)
    add_arg(args, "-v", v, 3)
    add_arg(args, "--gpu", gpu, False)
    add_arg(args, "--gpu-server", gpu_server, False)
    add_arg(args, "--mpi-runner", mpi_runner, "")
    add_arg(args, "--force-reuse", force_reuse, False)
    add_arg(args, "--remove-tmp-files", remove_tmp_files, True)

    # Expert
    add_arg(args, "--filter-hits", filter_hits, False)
    add_arg(args, "--sort-results", sort_results, 0)
    add_arg(args, "--create-lookup", create_lookup, False)
    add_arg(args, "--chain-alignments", chain_alignments, False)
    add_arg(args, "--merge-query", merge_query, True)
    add_arg(args, "--strand", strand, 1)
    add_arg(args, "--db-output", db_output, False)
    add_arg(args, "--write-lookup", write_lookup, False)

    # Execute command
    mmseqs_output = run_mmseqs_command(args)
    print(mmseqs_output.stdout)
    if mmseqs_output.stderr:
        print(mmseqs_output.stderr)
    print(f"MMseqs2 easy-search completed. Results saved to: {alignment_file_path}")
