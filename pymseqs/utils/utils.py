# pymseqs/utils/utils.py

import os
import inspect
from pathlib import Path
from typing import Any, Tuple, List

def get_caller_dir() -> str:
    """
    Get the directory of the script that's using this function.
    
    Returns:
        str: Absolute path to the directory containing the calling script
    """
    frame = inspect.currentframe()
    try:
        # Get the first frame (one level up)
        caller_frame = frame.f_back
        if caller_frame is None:
            return os.getcwd()
            
        # Check if we need to go up one more level
        caller_name = caller_frame.f_code.co_name
        if caller_name != '<module>':  # If we're in a function, go up one more level
            caller_frame = caller_frame.f_back
            
        # Get the full path of the calling script
        caller_file = caller_frame.f_code.co_filename
        # Return the directory containing the script
        return os.path.dirname(os.path.abspath(caller_file))
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

def add_twin_arg(
    args: List,
    flag: str,
    value: Tuple,
    defaults: Tuple,
    sep: str
):
    if value is not None and value != defaults:
        args.extend([flag, f"{value[0]}{sep}{value[1]}"])