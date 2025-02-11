from pathlib import Path
from typing import Union, List
from ..config.base import BaseConfig
from ..defaults import loader

DEFAULTS = loader.load("createdb")

class CreateDBConfig(BaseConfig):
    """
    Create a MMseqs2 database from a FASTA file and save it to the specified path prefix

    Parameters
    ----------
    `input_fasta` : Union[List[Union[str, Path]], Union[str, Path]]
    Path(s) to the input FASTA file(s). This can be:
        - A single string or Path object (e.g., `"input.fasta"` or `Path("input.fasta")`)
        - A list of strings or Path objects (e.g., `["input1.fasta", "input2.fasta"]` or 
          `[Path("input1.fasta"), Path("input2.fasta")]`)
        - A mixed list of strings and Path objects (e.g., `["input1.fasta", Path("input2.fasta")]`)

        If the paths are relative, they are resolved relative to the directory of the calling script
        If the paths are absolute, they are used as-is.

    `db_name` : Union[str, Path]
        Database path prefix, including the desired directory structure (e.g., `"output/dbs/mydb"`)
        - If the path is relative, it is resolved relative to the directory
        of the calling script
        - If the path is absolute, it is used as-is
        All necessary parent directories will be created automatically

    `dbtype` : int, optional
        Database type
        - 0: Auto-detect (default)
        - 1: Amino acid sequences
        - 2: Nucleotide sequences

    `shuffle` : bool, optional
        Shuffle the input database entries
        - True (default)
        - False

    `createdb_mode` : int, optional
        Database creation mode
        - 0: Copy data (default)
        - 1: Soft-link data and write a new index (only works with single-line FASTA/Q)

    `id_offset` : int, optional
        Numeric ID offset in the index file
        - 0 (default)

    `compressed` : bool, optional
        Compress the output files
        - True
        - False (default)

    `v` : int, optional
        Verbosity level of the output
        - 0: Quiet
        - 1: Errors only
        - 2: Errors and warnings
        - 3: Errors, warnings, and info (default)

    `write_lookup` : bool, optional
        Create a `.lookup` file mapping internal IDs to FASTA IDs
        - True (default)
        - False
    """

    def __init__(
        self,
        input_fasta: Union[List[Union[str, Path]], Union[str, Path]],
        db_name: Union[str, Path],
        dbtype: int = DEFAULTS["dbtype"]["default"],
        shuffle: bool = DEFAULTS["shuffle"]["default"],
        createdb_mode: int = DEFAULTS["createdb_mode"]["default"],
        id_offset: int = DEFAULTS["id_offset"]["default"],
        compressed: bool = DEFAULTS["compressed"]["default"],
        v: int = DEFAULTS["v"]["default"],
        write_lookup: bool = DEFAULTS["write_lookup"]["default"]
    ):
        self.input_fasta = input_fasta if isinstance(input_fasta, list) else [input_fasta]
        self.input_fasta = [Path(f) for f in self.input_fasta]
        self.db_name = Path(db_name)
        self.dbtype = dbtype
        self.shuffle = shuffle
        self.createdb_mode = createdb_mode
        self.id_offset = id_offset
        self.compressed = compressed
        self.v = v
        self.write_lookup = write_lookup

        self._defaults = DEFAULTS
        self._path_params = ['input_fasta', 'db_name']
        self._required_files = ['input_fasta']

    def validate(self) -> None:
        """
        Validate the configuration parameters.
        Raises a ValueError if any parameter is invalid.
        """
        if self.dbtype not in self._defaults["dbtype"]["choices"]:
            raise ValueError(f"dbtype is {self.dbtype} but must be one of {self._defaults['dbtype']['choices']}")
        
        if self.createdb_mode not in self._defaults["createdb_mode"]["choices"]:
            raise ValueError(f"createdb_mode is {self.createdb_mode} but must be one of {self._defaults['createdb_mode']['choices']}")
            
        if self.id_offset < 0:
            raise ValueError(f"id_offset is {self.id_offset} but must be non-negative")
            
        if self.v not in self._defaults["v"]["choices"]:
            raise ValueError(f"verbosity (v) is {self.v} but must be one of {self._defaults['v']['choices']}")
        
        if self.shuffle not in self._defaults["shuffle"]["choices"]:
            raise ValueError(f"shuffle is {self.shuffle} but must be one of {self._defaults['shuffle']['choices']}")
        
        if self.compressed not in self._defaults["compressed"]["choices"]:
            raise ValueError(f"compressed is {self.compressed} but must be one of {self._defaults['compressed']['choices']}")
        
        if self.write_lookup not in self._defaults["write_lookup"]["choices"]:
            raise ValueError(f"write_lookup is {self.write_lookup} but must be one of {self._defaults['write_lookup']['choices']}")
        
        # Optionally, re-check each FASTA file exists
        for fasta in self.input_fasta:
            if not Path(fasta).exists():
                raise ValueError(f"Input FASTA file does not exist: {fasta}")