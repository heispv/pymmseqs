import unittest
import subprocess
from pathlib import Path
import tempfile

from pymseqs.commands.createdb import createdb

class TestCreateDB(unittest.TestCase):
    def test_createdb_output_matches_cli(self):
        """
        Test that our createdb function produces identical output to mmseqs CLI
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            
            # Create test FASTA file
            fasta_content = ">seq1\nAAAA\n>seq2\nCCCC\n"
            input_fasta = tmp_path / "input_test.fasta"
            input_fasta.write_text(fasta_content)
            
            # Define output paths
            func_output = tmp_path / "func_output" / "mydb"
            cli_output = tmp_path / "cli_output" / "mydb"
            
            # 1. Run our Python function implementation
            createdb(
                input_fasta=input_fasta,
                db_name=func_output,
                write_lookup=1  # Match CLI default behavior
            )
            
            # 2. Run mmseqs CLI command
            # Create parent directory for CLI output
            cli_output.parent.mkdir(parents=True, exist_ok=True)
            
            subprocess.run(
                [
                    "mmseqs",
                    "createdb",
                    str(input_fasta),
                    str(cli_output),
                    "--write-lookup", "1"  # Match our function's default
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            
            # 3. Compare outputs
            # Get all generated files from both implementations
            func_files = set(func_output.parent.glob("mydb*"))
            cli_files = set(cli_output.parent.glob("mydb*"))
            
            # Check same number of files created
            self.assertEqual(
                len(func_files),
                len(cli_files),
                "Different number of output files generated"
            )
            
            # Compare each pair of files
            for func_file in func_files:
                filename = func_file.name
                cli_file = cli_output.parent / filename
                
                # Verify file exists in CLI output
                self.assertTrue(
                    cli_file.exists(),
                    f"File {filename} missing in CLI output"
                )
                
                # Compare file contents
                with open(func_file, "rb") as f1, open(cli_file, "rb") as f2:
                    self.assertEqual(
                        f1.read(),
                        f2.read(),
                        f"Content mismatch in {filename}"
                    )

if __name__ == "__main__":
    unittest.main()