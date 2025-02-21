import os
import subprocess
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
from setuptools.dist import Distribution

class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True

class build_py(_build_py):
    def run(self):
        # Ensure the build_lib/pymmseqs/bin directory exists
        target_dir = os.path.join(self.build_lib, 'pymmseqs', 'bin')
        self.mkpath(target_dir)
        # Download the binary to target_dir if it doesn’t already exist
        if not os.path.exists(os.path.join(target_dir, 'mmseqs')):
            self.download_mmseqs(target_dir)
        # Run the standard build_py process
        super().run()

    def download_mmseqs(self, target_dir):
        # Call the shell script with the target directory
        try:
            subprocess.check_call(['bash', 'scripts/download_mmseqs.sh', target_dir])
        except subprocess.CalledProcessError as e:
            raise RuntimeError("Failed to download MMseqs binary") from e

setup(
    name="pymmseqs",
    packages=find_packages(),
    package_data={
        "pymmseqs": ["bin/*"]
    },
    include_package_data=True,
    distclass=BinaryDistribution,
    cmdclass={'build_py': build_py},
)
