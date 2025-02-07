# pymseqs/commands/createdb.py

from pathlib import Path
from typing import Union, List

from pymseqs import run_mmseqs_command
from pymseqs.utils import (
    get_caller_dir,
    resolve_path,
    add_arg
)

def createdb(
    input_fasta: Union[List[Union[str, Path]], Union[str, Path]],
    db_name: Union[str, Path], 
    dbtype: int = 0,
    shuffle: bool = True,
    createdb_mode: int = 0, 
    id_offset: int = 0,
    compressed: bool = False,
    v: int = 3, 
    write_lookup: bool = True
) -> None:
    """
    Create a MMseqs2 database from a FASTA file and save it to the specified path prefix

    Parameters
    ----------
    **input_fasta** : Union[List[Union[str, Path]], Union[str, Path]]
    Path(s) to the input FASTA file(s). This can be:
        - A single string or Path object (e.g., `"input.fasta"` or `Path("input.fasta")`)
        - A list of strings or Path objects (e.g., `["input1.fasta", "input2.fasta"]` or 
          `[Path("input1.fasta"), Path("input2.fasta")]`)
        - A mixed list of strings and Path objects (e.g., `["input1.fasta", Path("input2.fasta")]`)

        If the paths are relative, they are resolved relative to the directory of the calling script
        If the paths are absolute, they are used as-is.

    **db_name** : Union[str, Path]
        Database path prefix, including the desired directory structure (e.g., `"output/dbs/mydb"`)
        - If the path is relative, it is resolved relative to the directory
        of the calling script
        - If the path is absolute, it is used as-is
        All necessary parent directories will be created automatically

    **dbtype** : int, optional
        Database type
        - 0: Auto-detect (default)
        - 1: Amino acid sequences
        - 2: Nucleotide sequences

    **shuffle** : bool, optional
        Shuffle the input database entries
        - True (default)
        - False

    **createdb_mode** : int, optional
        Database creation mode
        - 0: Copy data (default)
        - 1: Soft-link data and write a new index (only works with single-line FASTA/Q)

    **id_offset** : int, optional
        Numeric ID offset in the index file
        - 0 (default)

    **compressed** : bool, optional
        Compress the output files
        - True
        - False (default)

    **v** : int, optional
        Verbosity level of the output
        - 0: Quiet
        - 1: Errors only
        - 2: Errors and warnings
        - 3: Errors, warnings, and info (default)

    **write_lookup** : bool, optional
        Create a `.lookup` file mapping internal IDs to FASTA IDs
        - True (default)
        - False

    Returns
    -------
    None
        This function does not return any value, but creates files in the specified output directory.
    """

    # Get the directory of the calling script
    caller_dir = get_caller_dir()

    # Convert input paths to Path objects
    if not isinstance(input_fasta, list):
        input_fasta = [input_fasta]
    # Process all input files
    input_paths = []
    for file in input_fasta:
        # Convert to Path object and resolve path
        file_path = Path(file)
        resolved_path = resolve_path(file_path, caller_dir)
        # Check existence for each file
        if not resolved_path.exists():
            raise FileNotFoundError(f"Input file not found: {resolved_path}")
        input_paths.append(str(resolved_path))

    db_name_path = Path(db_name)
    db_name_path = resolve_path(db_name_path, caller_dir)
    
    args = [
        "createdb",
        *input_paths,
        str(db_name_path)
    ]
    
    # Define options with their corresponding current values and default values
    add_arg(args, "--dbtype", dbtype, 0)
    add_arg(args, "--shuffle", shuffle, True)
    add_arg(args, "--createdb-mode", createdb_mode, 0)
    add_arg(args, "--id-offset", id_offset, 0)
    add_arg(args, "--compressed", compressed, False)
    add_arg(args, "-v", v, 3)
    add_arg(args, "--write-lookup", write_lookup, True)
    
    mmseqs_output = run_mmseqs_command(args)
    
    if mmseqs_output.returncode == 0:
        print(mmseqs_output.stdout) if mmseqs_output.stdout else None
        print(f"Database path: {db_name_path}")
