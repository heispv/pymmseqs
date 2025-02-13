# pymmseqs/config/createdb_config.py

from pathlib import Path
from typing import Union, List

from .base import BaseConfig
from ..defaults import loader   
from ..utils import (
    get_caller_dir,
    add_arg,
    run_mmseqs_command
)

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
        dbtype: int = 0,
        shuffle: bool = True,
        createdb_mode: int = 0,
        id_offset: int = 0,
        compressed: bool = False,
        v: int = 3,
        write_lookup: bool = True
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
        self._path_params = [param for param, info in DEFAULTS.items() if info['type'] == 'path']
        self._required_files = [param for param, info in DEFAULTS.items()
                                if info['required'] and info['should_exist']]

    def validate(self) -> None:
        """
        Validate the configuration parameters.
        Raises a ValueError if any parameter is invalid.
        """
        for param_name, param_info in self._defaults.items():
            value = getattr(self, param_name)
            
            # Skip optional parameters with default values
            if not param_info['required'] and value == param_info['default']:
                continue
                
            # Validate choices if they exist
            if param_info['choices'] is not None:
                if value not in param_info['choices']:
                    raise ValueError(
                        f"{param_name} is {value} but must be one of {param_info['choices']}"
                    )
            
            # Validate numeric constraints
            if param_info['type'] == 'int' and param_name == 'id_offset' and value < 0:
                raise ValueError(f"id_offset is {value} but must be non-negative")
            
            # Validate required files exist
            if param_info['required'] and param_info['should_exist']:
                if isinstance(value, list):
                    for file_path in value:
                        if not Path(file_path).exists():
                            raise ValueError(f"Input file does not exist: {file_path}")
                elif not Path(value).exists():
                    raise ValueError(f"Input file does not exist: {value}")

    def run(self) -> None:
        # Get the directory of the calling script
        caller_dir = Path(get_caller_dir())

        # Use the config method to resolve all paths
        self.resolve_all_path(caller_dir)

            # Validate that all required files exist
        self.validate_required_files()
            
        # Create the command arguments
        args = [
            "createdb",
            *self.input_fasta,
            str(self.db_name)
        ]
        
        # Loop through all the optional parameters and add the arguments
        for param_name, param_info in self._defaults.items():
            if param_info['required']:
                continue
            
            cmd_param = f"--{param_name.replace('_', '-')}"
            
            current_value = getattr(self, param_name)
            default_value = param_info['default']
            
            add_arg(args, cmd_param, current_value, default_value)

        # Run the command
        mmseqs_output = run_mmseqs_command(args)
        
        if mmseqs_output.returncode == 0:
            if mmseqs_output.stdout:
                print(mmseqs_output.stdout)
            print(f"Database path: {self.db_name}")