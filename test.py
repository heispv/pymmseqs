from pymmseqs.wrapper import createdb
import os

os.makedirs('testing_folder', exist_ok=True)

createdb('input.fasta', 'testing_folder/db')