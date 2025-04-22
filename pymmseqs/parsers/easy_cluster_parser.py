# pymmseqs/parsers/easy_cluster_parser.py

from typing import Generator
import pandas as pd

from ..tools.easy_cluster_tools import (
    parse_fasta_clusters
)
from ..config import EasyClusterConfig

class EasyClusterParser:
    """
    A class for parsing the output of the EasyClusterConfig.
    """
    def __init__(
        self,
        config: EasyClusterConfig,
        seq_id_separator: str = "|",
        seq_id_index: int = 1
    ):
        """
        Parameters
        ----------
        config: EasyClusterConfig
            The configuration object for the EasyCluster command.
        seq_id_separator: str, optional
            The separator used in the FASTA headers to separate the sequence ID from other information.
            Default is "|".
        seq_id_index: int, optional
            The index of the sequence ID in the FASTA header.
            Default is 1.
        Note: It tries to extract the seq_id from the header using the separator and index, if it fails, it doesn't add the seq_id to the member.
        """
        self.cluster_prefix = config.cluster_prefix
        self.seq_id_separator = seq_id_separator
        self.seq_id_index = seq_id_index
    
    def to_list(self) -> list:
        """
        Parses a FASTA file containing clustered sequences and returns a list of dictionaries,
        where each dictionary represents a cluster.

        Returns:
        --------
        list of dict
            A list of dictionaries where each dictionary represents a single cluster with the following keys:
            - "rep": The representative sequence ID.
            - "members": List of member dictionaries in the cluster with the following keys:
                - "seq_id": Unique sequence identifier extracted from the header.
                    - If the header has format like ">seq_id|header", the seq_id is extracted from the header.
                - "header": Full FASTA header for the sequence.
                - "sequence": Nucleotide or protein sequence.

        When to Use:
        ------------
        - When you need to preserve the order of clusters as they appear in the file.
        - When you need to process all clusters at once and memory usage is not a concern.
        """
        return [
            {
            "rep": rep,
            "members": members
            }
            for rep, members in parse_fasta_clusters(f"{self.cluster_prefix}_all_seqs.fasta", self.seq_id_separator, self.seq_id_index)
        ]
    
    def to_pandas(self) -> pd.DataFrame:
        clusters = self.to_list()
        rows = []
        for cluster in clusters:
            rep = cluster["rep"]
            for member in cluster["members"]:
                # If the member has a seq_id, add it to the rows
                if "seq_id" in member:
                    rows.append({
                        "rep": rep,
                        "seq_id": member["seq_id"],
                        "header": member["header"],
                        "sequence": member["sequence"]
                    })
                else:
                    rows.append({
                        "rep": rep,
                        "header": member["header"],
                        "sequence": member["sequence"]
                    })
        return pd.DataFrame(rows).set_index('rep')
    
    def to_gen(self) -> Generator:
        """
        Generator that yields clusters one at a time from a FASTA file as dictionaries.

        Yields:
        -------
        dict
            A dictionary which represents a single cluster with the following keys:
            - "rep": The representative sequence ID.
            - "members": List of member dictionaries in the cluster with the following keys:
                - "seq_id": Unique sequence identifier extracted from the header.
                - "header": Full FASTA header for the sequence.
                - "sequence": Nucleotide or protein sequence.

        When to Use:
        ------------
        - When processing very large files where loading all clusters at once would consume too much memory.
        - When implementing streaming pipelines that process one cluster at a time.
        - When you need a dictionary format but want to avoid loading the entire dataset into memory.
        """
        for rep, members in parse_fasta_clusters(f"{self.cluster_prefix}_all_seqs.fasta", self.seq_id_separator, self.seq_id_index):
            yield {
                "rep": rep,
                "members": members
            }
    
    def to_path(self) -> list[str]:
        """
        Returns a list of file paths for the output files.

        Returns:
        --------
        list of str
        """
        return [
            f"{self.cluster_prefix}_all_seqs.fasta",
            f"{self.cluster_prefix}_cluster.tsv",
            f"{self.cluster_prefix}_rep_seqs.fasta",
        ]
