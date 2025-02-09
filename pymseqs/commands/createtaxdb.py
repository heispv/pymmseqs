# pymseqs/commands/createtaxdb.py

from pathlib import Path
from typing import Union, Optional

from pymseqs import run_mmseqs_command
from pymseqs.utils import (
    get_caller_dir,
    resolve_path,
    add_arg
)

def createtaxdb(
    sequence_db: Union[str, Path],
    tmp_dir: Union[str, Path],

    ncbi_tax_dump: Optional[Union[str, Path]] = None,
    tax_mapping_file: Optional[Union[str, Path]] = None,
    tax_mapping_mode: int = 0,
    tax_db_mode: int = 1,
    threads: int = 14,
    v: int = 3
) -> None:
    """
    Create a taxonomy database for a MMseqs2 sequence database

    Parameters
    ----------
    sequence_db : Union[str, Path]
        Path to the input MMseqs2 sequence database
        - If the path is relative, it is resolved relative to the directory
        of the calling script
        - If the path is absolute, it is used as-is

    tmp_dir : Union[str, Path]
        Path to the temporary directory
        - If the path is relative, it is resolved relative to the directory
        of the calling script
        - If the path is absolute, it is used as-is
        All necessary parent directories will be created automatically

    ncbi_tax_dump : Optional[Union[str, Path]], optional
        Path to the NCBI taxonomy dump directory
        The tax dump can be downloaded from:
        `ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`

    tax_mapping_file : Optional[Union[str, Path]], optional
        Path to the file mapping sequence identifiers to taxonomical identifiers

    tax_mapping_mode : int, optional
        Mode for mapping taxonomy based on sequence database
        - 0: .lookup file (default)
        - 1: .source file

    tax_db_mode : int, optional
        Mode for creating taxonomy database
        - 0: .dmp flat files (human readable)
        - 1: binary dump (faster reading) (default)

    threads : int, optional
        Number of CPU cores to use
        - Default: 14

    v : int, optional
        Verbosity level of the output
        - 0: quiet
        - 1: + errors
        - 2: + warnings
        - 3: + info (default)

    Returns
    -------
    None
        This function does not return any value, but creates taxonomy database files
        in the same directory as the sequence database.
    """
    # Get the directory of the calling script
    caller_dir = get_caller_dir()

    # Convert and resolve paths
    sequence_db_path = resolve_path(Path(sequence_db), caller_dir)
    tmp_dir_path = resolve_path(Path(tmp_dir), caller_dir)

    # Check if sequence database exists
    if not sequence_db_path.exists():
        raise FileNotFoundError(f"Sequence database not found: {sequence_db_path}")

    # Create tmp directory if it doesn't exist
    tmp_dir_path.mkdir(parents=True, exist_ok=True)

    args = [
        "createtaxdb",
        str(sequence_db_path),
        str(tmp_dir_path)
    ]

    # Handle optional NCBI tax dump directory
    if ncbi_tax_dump is not None:
        ncbi_tax_dump_path = resolve_path(Path(ncbi_tax_dump), caller_dir)
        if not ncbi_tax_dump_path.exists():
            raise FileNotFoundError(f"NCBI tax dump directory not found: {ncbi_tax_dump_path}")
        add_arg(args, "--ncbi-tax-dump", str(ncbi_tax_dump_path))

    # Handle optional taxonomy mapping file
    if tax_mapping_file is not None:
        tax_mapping_path = resolve_path(Path(tax_mapping_file), caller_dir)
        if not tax_mapping_path.exists():
            raise FileNotFoundError(f"Taxonomy mapping file not found: {tax_mapping_path}")
        add_arg(args, "--tax-mapping-file", str(tax_mapping_path))

    # Add other arguments
    add_arg(args, "--tax-mapping-mode", tax_mapping_mode, 0)
    add_arg(args, "--tax-db-mode", tax_db_mode, 1)
    add_arg(args, "--threads", threads, 14)
    add_arg(args, "-v", v, 3)

    mmseqs_output = run_mmseqs_command(args)

    if mmseqs_output.returncode == 0:
        print(mmseqs_output.stdout) if mmseqs_output.stdout else None
        print(f"Taxonomy database created for: {sequence_db_path}")
