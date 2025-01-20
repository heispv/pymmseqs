# pymseqs/commands/search.py

from pymseqs import run_mmseqs_command

def search(query, db, result, threads=None, **kwargs):
    """
    Perform a search using MMseqs2.
    
    Parameters:
        query (str): Path to the query FASTA file.
        db (str): Path to the MMseqs2 database.
        result (str): Path to the output result file.
        threads (int, optional): Number of threads to use.
        **kwargs: Additional MMseqs2 search parameters.
    
    Returns:
        stdout (str): Output from the mmseqs2 command.
    """
    args = ['search', query, db, result]
    
    # Example mapping for threads; extend this as needed
    if threads is not None:
        args.extend(['--threads', str(threads)])
    
    # Handle additional search options if any
    for key, value in kwargs.items():
        args.extend([f'--{key}', str(value)])
    
    return run_mmseqs_command(args)