# pymseqs/utils/utils.py

import os
import inspect

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