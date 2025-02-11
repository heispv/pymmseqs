# pymseqs/commands/createdb.py

from pathlib import Path

from ..runner import run_mmseqs_command
from ..config import CreateDBConfig
from ..utils import (
    get_caller_dir,
    add_arg
)

def createdb(config: CreateDBConfig) -> None:
    """
    Create a MMseqs2 database from a FASTA file and save it to the specified path prefix

    Parameters
    ----------
    config : CreateDBConfig
        Configuration object containing all parameters for the createdb command.

    Returns
    -------
    None
        This function does not return any value, but creates files in the specified output directory.
    """

    # Get the directory of the calling script
    caller_dir = Path(get_caller_dir())

    # Use the config method to resolve all paths
    config.resolve_all_path(caller_dir)

    # Validate that all required files exist
    config.validate_required_files()
    
    # Create the command arguments
    args = [
        "createdb",
        *config.input_fasta,
        str(config.db_name)
    ]
    
    # Loop through all the optional parameters and add the arguments
    for param_name, param_info in config._defaults.items():
        if not param_info['optional']:
            continue

        cmd_param = f"--{param_name.replace('_', '-')}"
        
        current_value = getattr(config, param_name)
        default_value = param_info['default']
        
        add_arg(args, cmd_param, current_value, default_value)

    # Run the command
    mmseqs_output = run_mmseqs_command(args)
    
    if mmseqs_output.returncode == 0:
        if mmseqs_output.stdout:
            print(mmseqs_output.stdout)
        print(f"Database path: {config.db_name}")
