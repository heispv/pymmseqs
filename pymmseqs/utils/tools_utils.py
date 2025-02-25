# pymmseqs/utils/tools_utils.py
import pandas as pd
from csv import Sniffer

def to_superscript(exp):
    superscripts = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
    return str(exp).translate(superscripts)


def has_header(file_path):
    """
    Determines if a CSV or TSV file has a header row.
    """
    try:
        with open(file_path, 'r', newline='') as file:
            sample = file.read(1024)
            if not sample.strip():
                return False
            
            sniffer = Sniffer()
            dialect = sniffer.sniff(sample)
            file.seek(0)
            reader = pd.read_csv(file, sep=dialect.delimiter, nrows=3)
            if reader.shape[0] == 0:
                return False
            
            headers = reader.columns
            for header in headers:
                if 'evalue' in header.lower().replace('-', ''):
                    return True
            
            row = reader.iloc[0]
            return any('e' in str(cell).lower() for cell in row)
    except Exception:
        return False
