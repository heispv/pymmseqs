# pymseqs/commands/createdb.py

from pathlib import Path
from typing import Union
import os

from pymseqs import run_mmseqs_command
from pymseqs.utils import get_caller_dir

def createdb(
    input_fasta: Union[str, Path],
    db_name: Union[str, Path], 
    dbtype: int = 0,
    shuffle: int = 1,
    createdb_mode: int = 0, 
    id_offset: int = 0,
    compressed: int = 0,
    verbosity: int = 3, 
    write_lookup: int = 1
) -> None:
    """
    Create a MMseqs2 database from a FASTA file and save it to the specified path prefix.

    Parameters
    ----------
    input_fasta : Union[str, Path]
        Path to the input FASTA file. If the path is relative, it is resolved relative to
        the directory of the calling script. If the path is absolute, it is used as-is.

    db_name : Union[str, Path]
        Database path prefix, including the desired directory structure (e.g., 
        "output/dbs/mydb").
        - If the path is relative, it is resolved relative to the directory
        of the calling script.
        - If the path is absolute, it is used as-is.
        All necessary parent directories will be created automatically.

    dbtype : int, optional
        Database type. Options:
        - 0: Auto-detect (default)
        - 1: Amino acid sequences
        - 2: Nucleotide sequences

    shuffle : int, optional
        Whether to shuffle the input database entries. Options:
        - 0: Disabled
        - 1: Enabled (default)

    createdb_mode : int, optional
        Database creation mode. Options:
        - 0: Copy data (default)
        - 1: Soft-link data and write a new index (only works with single-line FASTA/Q)

    id_offset : int, optional
        Numeric ID offset in the index file. Default is 0.

    compressed : int, optional
        Whether to compress the output files. Options:
        - 0: Uncompressed (default)
        - 1: Compressed

    verbosity : int, optional
        Verbosity level of the output. Options:
        - 0: Quiet
        - 1: Errors only
        - 2: Errors and warnings
        - 3: Errors, warnings, and info (default)

    write_lookup : int, optional
        Whether to create a `.lookup` file mapping internal IDs to FASTA IDs. Options:
        - 0: Disabled
        - 1: Enabled (default)

    Returns
    -------
    None
        This function does not return any value, but creates files in the specified output directory.
    """

    # Get the directory of the calling script
    caller_dir = get_caller_dir()

    # Convert input paths to Path objects
    input_fasta_path = Path(input_fasta)
    db_name_path = Path(db_name)
    
    # Resolve input_fasta relative to caller_dir if not absolute
    if not input_fasta_path.is_absolute():
        input_fasta_path = caller_dir / input_fasta_path
    
    # Resolve db_name relative to caller_dir if not absolute
    if not db_name_path.is_absolute():
        db_name_path = caller_dir / db_name_path
    
    # Create parent directory for db_name
    db_parent = db_name_path.parent
    os.makedirs(db_parent, exist_ok=True)
    
    # Ensure input file exists
    if not input_fasta_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_fasta_path}")
    
    args = ['createdb', str(input_fasta_path), str(db_name_path)]
    
    # Define options with their corresponding current values and default values
    options = [
        ('--dbtype', dbtype, 0),
        ('--shuffle', int(shuffle), 1),
        ('--createdb-mode', createdb_mode, 0),
        ('--id-offset', id_offset, 0),
        ('--compressed', compressed, 0),
        ('-v', verbosity, 3),
        ('--write-lookup', write_lookup, 1),
    ]
    
    # Append options only if they differ from their default values
    for option, value, default in options:
        if value != default:
            args.extend([option, str(value)])
    
    mmseqs_output = run_mmseqs_command(args)
    print(mmseqs_output)
