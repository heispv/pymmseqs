# setup.py
import os
import sys

# PEP 517 workaround: Add project root to Python path
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py

from setup_utils import download_mmseqs2

class CustomBuildCommand(build_py):
    """Custom build command to download MMseqs2 before building."""
    def run(self):
        install_dir = os.path.join(self.build_lib, "pymmseqs", "bin")
        download_mmseqs2(install_dir)
        super().run()

setup(
    packages=find_packages(),
    # package_data={"pymseqs": ["bin/*"],
    #               "pymseqs.defaults": ["*.yaml"]},
    include_package_data=True,
    cmdclass={"build_py": CustomBuildCommand},
)
