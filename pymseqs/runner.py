# pymseqs/runner.py

import subprocess
from pymseqs import get_mmseqs_binary

def run_mmseqs_command(args, capture_output=True):
    """
    Run an mmseqs2 command with the given arguments.
    Raises RuntimeError if the command fails.
    Returns the command's stdout.
    """
    binary = get_mmseqs_binary()
    cmd = [binary] + args
    print(f"Executing command: {' '.join(cmd)}")  # Debug statement
    
    result = subprocess.run(cmd, capture_output=capture_output, text=True)
    if result.returncode != 0:
        print(f"mmseqs2 stderr: {result.stderr}")  # Debug statement
        raise RuntimeError(f"mmseqs2 failed: {result.stderr}")
    
    return result.stdout
