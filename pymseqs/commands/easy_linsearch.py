# pymseqs/commands/easy_linsearch.py

from pathlib import Path
from typing import Union, List, Tuple

from pymseqs import run_mmseqs_command
from pymseqs.utils import (
    get_caller_dir,
    resolve_path,
    add_arg,
    add_twin_arg
)


def easy_linsearch(
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
    mask: bool = True,
    mask_prob: float = 0.9,
    mask_lower_case: bool = False,
    split_memory_limit: str = "0",
    
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

    # K-mer matcher parameters
    kmer_per_seq: int = 21,
    kmer_per_seq_scale: Tuple[str, str] = ("aa:0.0", "nucl:0.2"),
    pick_n_sim_kmer: int = 1,
    result_direction: int = 1,
    
    # Profile parameters
    pca: float = 0.0,
    pcb: float = 0.0,
    
    # Misc parameters
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
    search_type: int = 0,
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
    
    # Common parameters
    sub_mat: Tuple[str, str] = ("aa:blosum62.out", "nucl:nucleotide.out"),
    max_seq_len: int = 65535,
    db_load_mode: int = 0,
    threads: int = 14,
    compressed: bool = False,
    v: int = 3,
    mpi_runner: str = "",
    force_reuse: bool = False,
    remove_tmp_files: bool = True,
    
    # Expert parameters
    create_lookup: bool = False,
    chain_alignments: bool = False,
    merge_query: bool = True,
    db_output: bool = False,
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

    **split_memory_limit** : str, optional
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
    
    K-mer Matcher Parameters
    ------------------
    **kmer_per_seq** : int, optional
        Number of k-mers per sequence.
        - 21 (default)

    **kmer_per_seq_scale** : Tuple[str, str], optional
        Scale k-mer per sequence based on sequence length as (kmer-pers-seq val + scale * seq-len)
        - ("aa:0.0", "nucl:0.2") (default)
    
    **pick_n_sim_kmer** : int, optional
        Add N similar k-mers to search
        - True (default)
        - False
    
    **result_direction** : int, optional
        Results is
        - 0: query
        - 1: target centric (default)

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

    Misc Parameters
    ---------------
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

    **search_type** : int, optional
        Search mode:
        - 0: auto (default)
        - 1: amino
        - 2: translated
        - 3: nucleotide
        - 4: translated alignment
    
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

    Common Parameters
    ----------------
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

    **v** : int, optional
        Output verbosity
        - 0: quiet
        - 1: +errors
        - 2: +warnings
        - 3: +info (default)

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

    **db_output** : bool, optional
        Return the result as DB instead of a text file
        - True: DB output
        - False: Text file output (default)

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
        threads=8
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

    # Validate input files
    for fasta_path, name in [(query_fasta_path, "Query"), (target_fasta_or_db_path, "Target")]:
        if not fasta_path.exists():
            raise FileNotFoundError(f"{name} FASTA file not found: {fasta_path}")

    # Build base command
    args = [
        "easy-linsearch",
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
    add_arg(args, "--split-memory-limit", split_memory_limit, "0")
    add_arg(args, "--mask", mask, True)
    add_arg(args, "--mask-prob", mask_prob, 0.9)
    add_arg(args, "--mask-lower-case", mask_lower_case, False)
    add_arg(args, "--disk-space-limit", split_memory_limit, "0")

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

    # K-mer Matcher
    add_arg(args, "--kmer-per-seq", kmer_per_seq, 21)
    add_twin_arg(args, "--kmer-per-seq-scale", kmer_per_seq_scale, ("aa:0.0", "nucl:0.2"), ",")
    add_arg(args, "--pick-n-sim-kmer", pick_n_sim_kmer, 1)
    add_arg(args, "--result-direction", result_direction, 1)

    # Profile
    add_arg(args, "--pca", pca, 0.0)
    add_arg(args, "--pcb", pcb, 0.0)

    # Misc
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
    add_arg(args, "--search-type", search_type, 0)
    add_arg(args, "--format-mode", format_mode, 0)
    add_arg(
    args, "--format-output", ",".join(map(str, format_output)),
    "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
    )
    
    # Common
    add_twin_arg(args, "--sub-mat", sub_mat, ("aa:blosum62.out", "nucl:nucleotide.out"), ",")
    add_arg(args, "--max-seq-len", max_seq_len, 65535)
    add_arg(args, "--db-load-mode", db_load_mode, 0)
    add_arg(args, "--threads", threads, 14)
    add_arg(args, "--compressed", compressed, False)
    add_arg(args, "-v", v, 3)
    add_arg(args, "--mpi-runner", mpi_runner, "")
    add_arg(args, "--force-reuse", force_reuse, False)
    add_arg(args, "--remove-tmp-files", remove_tmp_files, True)

    # Expert
    add_arg(args, "--create-lookup", create_lookup, False)
    add_arg(args, "--chain-alignments", chain_alignments, False)
    add_arg(args, "--merge-query", merge_query, True)
    add_arg(args, "--db-output", db_output, False)

    # Execute command
    mmseqs_output = run_mmseqs_command(args)
    print(mmseqs_output.stdout)
    if mmseqs_output.stderr:
        print(mmseqs_output.stderr)
    print(f"MMseqs2 easy-search completed. Results saved to: {alignment_file_path}")
