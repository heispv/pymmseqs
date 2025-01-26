# pymseqs/utils/utils.py

import os
import inspect
from pathlib import Path
from typing import Any, Tuple, List

def get_caller_dir() -> str:
    """
    Get the directory of the calling script.
    
    Returns:
        str: Absolute path to the directory containing the calling script
    """
    frame = inspect.currentframe()
    # Get the frame of the calling function
    caller_frame = frame.f_back.f_back  # Two levels up because of the additional function call
    if caller_frame is not None:
        caller_file = caller_frame.f_code.co_filename
        caller_dir = os.path.dirname(os.path.abspath(caller_file))
    else:
        caller_dir = os.getcwd()
    return caller_dir

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