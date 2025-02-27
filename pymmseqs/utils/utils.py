# pymmseqs/utils/utils.py

import os
import inspect
from pathlib import Path
from typing import Any, Tuple, List

import os
import sys
import inspect
from pathlib import Path
from IPython import get_ipython

def get_caller_dir() -> Path:
    """
    Get the directory of the script that's using this function.
    
    For .py files, traverses the call stack to find the first frame outside
    the pymmseqs package, presumed to be the user's code. For .ipynb files
    (Jupyter notebooks), returns the current working directory, which is
    typically the directory containing the notebook unless changed.
    
    Returns:
        Path: Absolute path to the directory containing the calling script
    """
    # Check if running in a Jupyter notebook
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            # Jupyter notebook detected; return current working directory
            return Path(os.getcwd())
    except NameError:
        # Not in an IPython environment; proceed with stack traversal
        pass
    
    # Get the full call stack
    frame = inspect.currentframe()
    try:
        # Get package path to identify frames within pymmseqs
        pymmseqs_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
        # Start from the immediate caller
        caller_frame = frame.f_back
        
        # Traverse up the stack until we find a frame outside pymmseqs
        while caller_frame:
            caller_file = caller_frame.f_code.co_filename
            
            # If the frame is not from within pymmseqs package or standard library
            if (not caller_file.startswith(pymmseqs_path) and 
                not caller_file.startswith(sys.prefix) and
                not caller_file == '<string>'):  # Ignore REPL or eval frames
                # Found a frame outside pymmseqs - likely the user's code
                return Path(os.path.dirname(os.path.abspath(caller_file)))
            
            # Move up to the next frame
            caller_frame = caller_frame.f_back
        
        # If no suitable frame is found, return current working directory
        return Path(os.getcwd())
    finally:
        # Clean up the frame to prevent memory leaks
        del frame

def resolve_path(
    path: Path,
    caller_dir: Path
) -> Path:
    """Resolves a path relative to `caller_dir` if not absolute and ensures its parent directory exists.

    Parameters
    ----------
    path : Path
        Input path (relative or absolute).
    caller_dir : Path
        Base directory for resolving relative paths.

    Returns
    -------
    Path
        Resolved absolute path. Parent directory is created if it doesn't exist.
    """
    path = Path(path)
    # Resolve relative path if not absolute
    if not path.is_absolute():
        path = caller_dir / path

    # Optionally create the parent directory
    os.makedirs(path.parent, exist_ok=True)

    return path

def add_arg(
    args: List,
    flag: str,
    value: Any,
    default: Any,
):
    if value != default:
        if isinstance(value, bool):
            args.extend([flag, "1" if value else "0"])
        else:
            args.extend([flag, str(value)])
