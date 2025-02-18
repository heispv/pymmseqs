from pathlib import Path
from typing import Union, List

from pymmseqs.config.base import BaseConfig
from pymmseqs.defaults import loader
from pymmseqs.utils import (
    get_caller_dir,
    run_mmseqs_command
)

DEFAULTS = loader.load("easy_cluster")

class EasyClusterConfig(BaseConfig):
    def __init__(
        self,
        # Required parameters
        fasta_files: Union[str, Path, List[Union[str, Path]]],
        cluster_prefix: Union[str, Path],
        tmp_dir: Union[str, Path],

        # Prefilter parameters
        seed_sub_mat: str = "aa:VTML80.out,nucl:nucleotide.out",
        s: float = 4.0,
        k: int = 0,
        target_search_mode: int = 0,
        k_score: str = "seq:2147483647,prof:2147483647",
        alph_size: str = "aa:21,nucl:5",
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
        mask_n_repeat: int = 0,
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
        gap_open: str = "aa:11,nucl:5",
        gap_extend: str = "aa:1,nucl:2",
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
        kmer_per_seq_scale: str = "aa:0.0,nucl:0.2",
        adjust_kmer_len: bool = False,
        hash_shift: int = 67,
        include_only_extendable: bool = False,
        ignore_multi_kmer: bool = False,

        # Profile parameters
        pca: float = None,
        pcb: float = None,

        # Misc parameters
        taxon_list: str = "",
        rescore_mode: int = 0,
        dbtype: int = 0,
        shuffle: bool = True,
        createdb_mode: int = 1,
        id_offset: int = 0,

        # Common parameters
        sub_mat: str = "aa:blosum62.out,nucl:nucleotide.out",
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
    ):
        """
        Configuration for MMseqs2 easy-cluster command.
        
        For detailed parameter descriptions, see the easy_cluster.yaml file.
        """
        # Required parameters
        self.fasta_files = fasta_files if isinstance(fasta_files, list) else [fasta_files]
        self.cluster_prefix = cluster_prefix
        self.tmp_dir = tmp_dir

        # Prefilter parameters
        self.seed_sub_mat = seed_sub_mat
        self.s = s
        self.k = k
        self.target_search_mode = target_search_mode
        self.k_score = k_score
        self.alph_size = alph_size
        self.max_seqs = max_seqs
        self.split = split
        self.split_mode = split_mode
        self.split_memory_limit = split_memory_limit
        self.comp_bias_corr = comp_bias_corr
        self.comp_bias_corr_scale = comp_bias_corr_scale
        self.diag_score = diag_score
        self.exact_kmer_matching = exact_kmer_matching
        self.mask = mask
        self.mask_prob = mask_prob
        self.mask_lower_case = mask_lower_case
        self.mask_n_repeat = mask_n_repeat
        self.min_ungapped_score = min_ungapped_score
        self.add_self_matches = add_self_matches
        self.spaced_kmer_mode = spaced_kmer_mode
        self.spaced_kmer_pattern = spaced_kmer_pattern
        self.local_tmp = local_tmp

        # Alignment parameters
        self.c = c
        self.cov_mode = cov_mode
        self.a = a
        self.alignment_mode = alignment_mode
        self.alignment_output_mode = alignment_output_mode
        self.wrapped_scoring = wrapped_scoring
        self.e = e
        self.min_seq_id = min_seq_id
        self.min_aln_len = min_aln_len
        self.seq_id_mode = seq_id_mode
        self.alt_ali = alt_ali
        self.max_rejected = max_rejected
        self.max_accept = max_accept
        self.score_bias = score_bias
        self.realign = realign
        self.realign_score_bias = realign_score_bias
        self.realign_max_seqs = realign_max_seqs
        self.corr_score_weight = corr_score_weight
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.zdrop = zdrop

        # Clustering parameters
        self.cluster_mode = cluster_mode
        self.max_iterations = max_iterations
        self.similarity_type = similarity_type
        self.single_step_clustering = single_step_clustering
        self.cluster_steps = cluster_steps
        self.cluster_reassign = cluster_reassign

        # K-mer matcher parameters
        self.weights = weights
        self.cluster_weight_threshold = cluster_weight_threshold
        self.kmer_per_seq = kmer_per_seq
        self.kmer_per_seq_scale = kmer_per_seq_scale
        self.adjust_kmer_len = adjust_kmer_len
        self.hash_shift = hash_shift
        self.include_only_extendable = include_only_extendable
        self.ignore_multi_kmer = ignore_multi_kmer

        # Profile parameters
        self.pca = pca
        self.pcb = pcb

        # Misc parameters
        self.taxon_list = taxon_list
        self.rescore_mode = rescore_mode
        self.dbtype = dbtype
        self.shuffle = shuffle
        self.createdb_mode = createdb_mode
        self.id_offset = id_offset

        # Common parameters
        self.sub_mat = sub_mat
        self.max_seq_len = max_seq_len
        self.db_load_mode = db_load_mode
        self.threads = threads
        self.compressed = compressed
        self.v = v
        self.remove_tmp_files = remove_tmp_files
        self.force_reuse = force_reuse
        self.mpi_runner = mpi_runner

        # Expert parameters
        self.filter_hits = filter_hits
        self.sort_results = sort_results
        self.write_lookup = write_lookup

        self._defaults = DEFAULTS
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']

    def _validate(self) -> None:
        self._check_required_files()
        self._validate_choices()

        # Validate numerical ranges
        if not (0 <= self.comp_bias_corr_scale <= 1):
            raise ValueError("comp_bias_corr_scale must be between 0 and 1")
        if not (1.0 <= self.s <= 7.5):
            raise ValueError("Sensitivity (-s) must be between 1.0 and 7.5")
        if not (0.0 <= self.min_seq_id <= 1.0):
            raise ValueError("min_seq_id must be between 0.0 and 1.0")
        if not (0.0 <= self.mask_prob <= 1.0):
            raise ValueError("mask_prob must be between 0.0 and 1.0")

    def run(self) -> None:
        """Execute the MMseqs2 easy-cluster command with the current configuration."""
        # Get the directory of the calling script
        caller_dir = get_caller_dir()

        # Resolve all paths
        self._resolve_all_path(caller_dir)

        # Validate the configuration
        self._validate()

        # Get command arguments and run the command
        args = self._get_command_args("easy-cluster")
        print(args)
        mmseqs_output = run_mmseqs_command(args)

        if mmseqs_output.returncode == 0:
            if mmseqs_output.stdout:
                print(mmseqs_output.stdout)
            if mmseqs_output.stderr:
                print(mmseqs_output.stderr)
            print(f"Clustering results saved to: {self.cluster_prefix}")
        else:
            raise RuntimeError(f"MMseqs2 clustering failed: {mmseqs_output.stderr}")