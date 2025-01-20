# pymseqs/commands/createdb.py

import inspect
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
    Create a MMseqs2 database from a FASTA file and save it to a Path.
    Paths are resolved relative to the calling script's directory.
    
    Parameters:
        input_fasta (str): Path to the input FASTA file.
        db_name (str): Desired name for the created database.
        dbtype (int, optional): Database type (0: auto, 1: amino acid, 2: nucleotides), default is [0].
        shuffle (bool, optional): Shuffle input database, (0: False, 1: True), default is [1].
        createdb_mode (int, optional): Createdb mode (0: copy data, 1: soft link data and write new index (works only with single line fasta/q)), default is [0].
        id_offset (int, optional): Numeric ids in index file are offset by this value, default is [0].
        compressed (int, optional): Write compressed output (0: no, 1: yes), default is [0].
        verbosity (int, optional): Verbosity level (0: quiet, 1: +errors, 2: +warnings, 3: +info), default is [3].
        write_lookup (int, optional): Write .lookup file containing mapping from internal id, fasta id and file number (0: no, 1: yes), default is [0].
    
    Returns:
        None
    """
    # Get the directory of the calling script
    caller_dir = get_caller_dir()
    output_dir = Path(caller_dir) / 'output'

    os.makedirs(output_dir, exist_ok=True)
    
    # Convert input paths to Path objects
    input_fasta_path = Path(input_fasta)
    db_name_path = Path(db_name)
    

    # If the paths are not absolute, make them relative to the caller's directory
    if not input_fasta_path.is_absolute():
        input_fasta_path = caller_dir / input_fasta_path
    if not db_name_path.is_absolute():
        db_name_path = output_dir / db_name_path
    
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
