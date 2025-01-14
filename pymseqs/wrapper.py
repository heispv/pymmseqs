# pymseqs/wrapper.py

import os
import subprocess
import platform
from pathlib import Path

def get_mmseqs_binary():
    """
    Retrieve the path to the mmseqs2 binary.
    Allows overriding via the MMSEQS2_PATH environment variable.
    """
    # Check if the MMSEQS2_PATH environment variable is set
    custom_path = os.getenv('MMSEQS2_PATH')
    if custom_path:
        if os.path.exists(custom_path):
            return custom_path
        else:
            raise FileNotFoundError(f"mmseqs2 binary specified by MMSEQS2_PATH does not exist: {custom_path}")
    
    # Default behavior: locate within the site-packages bin directory
    from sysconfig import get_path
    system = platform.system()
    binary_name = 'mmseqs.exe' if system == 'Windows' else 'mmseqs'
    binary_path = os.path.join(get_path('purelib'), 'pymseqs', 'bin', binary_name)
    
    if not os.path.exists(binary_path):
        raise FileNotFoundError(f"mmseqs2 binary not found at {binary_path}. Please ensure it is installed correctly.")
    
    return binary_path

def run_mmseqs_command(args, capture_output=True):
    """
    Run a mmseqs2 command with the given arguments.
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
    print(result.stdout)

def createdb(input_fasta, db_name, **kwargs):
    """
    Create a MMseqs2 database from a FASTA file.
    
    Parameters:
        input_fasta (str): Path to the input FASTA file.
        db_name (str): Desired name for the created database.
        **kwargs: Additional MMseqs2 parameters (e.g., threads).
    
    Returns:
        stdout (str): Output from the mmseqs2 command.
    """
    args = ['createdb', input_fasta, db_name]
    # Append additional arguments
    for key, value in kwargs.items():
        args.extend([f'--{key}', str(value)])
    return run_mmseqs_command(args)

def search(query, db, result, **kwargs):
    """
    Perform a search using MMseqs2.
    
    Parameters:
        query (str): Path to the query FASTA file.
        db (str): Path to the MMseqs2 database.
        result (str): Path to the output result file.
        **kwargs: Additional MMseqs2 parameters (e.g., threads).
    
    Returns:
        stdout (str): Output from the mmseqs2 command.
    """
    args = ['search', query, db, result]
    # Append additional arguments
    for key, value in kwargs.items():
        args.extend([f'--{key}', str(value)])
    return run_mmseqs_command(args)

def main():
    """
    Entry point for the pymseqs CLI.
    Parses command-line arguments and invokes the appropriate function.
    """
    import argparse

    parser = argparse.ArgumentParser(description='pymseqs CLI')
    subparsers = parser.add_subparsers(dest='command', required=True, help='Available commands')

    # Subparser for 'createdb' command
    createdb_parser = subparsers.add_parser('createdb', help='Create a MMseqs2 database')
    createdb_parser.add_argument('input_fasta', help='Path to the input FASTA file')
    createdb_parser.add_argument('db_name', help='Name for the created database')
    # Add optional arguments (e.g., threads)
    createdb_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')

    # Subparser for 'search' command
    search_parser = subparsers.add_parser('search', help='Perform a MMseqs2 search')
    search_parser.add_argument('query', help='Path to the query FASTA file')
    search_parser.add_argument('db', help='Path to the MMseqs2 database')
    search_parser.add_argument('result', help='Path to the output result file')
    # Add optional arguments (e.g., threads)
    search_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')

    args = parser.parse_args()

    if args.command == 'createdb':
        try:
            output = createdb(args.input_fasta, args.db_name, threads=args.threads)
            print(output)
        except Exception as e:
            print(f"Error during createdb: {e}")
    elif args.command == 'search':
        try:
            output = search(args.query, args.db, args.result, threads=args.threads)
            print(output)
        except Exception as e:
            print(f"Error during search: {e}")
    else:
        parser.print_help()
