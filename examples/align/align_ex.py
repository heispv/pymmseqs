from pymmseqs.config import (
    CreateDBConfig,
    SearchConfig,
    AlignConfig,
    ConvertAlisConfig
)

query_db = CreateDBConfig(
    fasta_file="human.fasta",
    sequence_db="output/query_db"
)

target_db = CreateDBConfig(
    fasta_file="human.fasta",
    sequence_db="output/target_db"
)

search_db = SearchConfig(
    query_db="output/query_db",
    target_db="output/target_db",
    alignment_db="output/result_db",
    tmp_dir="output/tmp",
    min_seq_id=0.95,
    max_seqs=1000000,

)

align_db = AlignConfig(
    query_db="output/query_db",
    target_db="output/target_db",
    result_db="output/result_db",
    alignment_db="output/alignment_db",
    alignment_mode=3
)

convert_alis = ConvertAlisConfig(
    query_db="output/query_db",
    target_db="output/target_db",
    alignment_db="output/alignment_db",
    alignment_file="output/alignment.m8",
    format_mode=4
)

query_db.run()
target_db.run()
search_db.run()
align_db.run()
convert_alis.run()
