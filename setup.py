import os
import sys
import io
import platform
import shutil
import tarfile
import zipfile
from urllib.error import URLError
from urllib.request import urlopen
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py

MMSEQS_VERSION = "16-747c6"

def get_mmseqs_download_info():
    """Determine the appropriate MMseqs2 binary URL and filename for the current platform."""
    system = platform.system()
    machine = platform.machine().lower()

    base_url = "https://github.com/soedinglab/mmseqs2/releases/download"

    if system == "Linux":
        if machine in ("x86_64", "amd64"):
            return f"{base_url}/{MMSEQS_VERSION}/mmseqs-linux-avx2.tar.gz", "mmseqs"
        elif machine == "aarch64":
            return f"{base_url}/{MMSEQS_VERSION}/mmseqs-linux-arm64.tar.gz", "mmseqs"
    elif system == "Darwin":
        return f"{base_url}/{MMSEQS_VERSION}/mmseqs-osx-universal.tar.gz", "mmseqs"
    elif system == "Windows":
        return f"{base_url}/{MMSEQS_VERSION}/mmseqs-win64.zip", "mmseqs.exe"

def download_mmseqs_binary(url):
    """Download MMseqs2 binary from the given URL."""
    print(f"Downloading MMseqs2 from {url}")
    try:
        with urlopen(url) as response:
            return response.read()
    except URLError as e:
        raise RuntimeError(f"Failed to download MMseqs2: {e}") from e

def extract_mmseqs_binary(data, url, install_dir, binary_name):
    """Extract MMseqs2 binary from downloaded data to install directory."""
    try:
        if url.endswith(".tar.gz"):
            with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tar:
                member_path = "mmseqs/bin/mmseqs"  # Adjusted path for Linux/macOS
                member = tar.getmember(member_path)
                member.name = os.path.basename(member.name)
                tar.extract(member, path=install_dir)
        elif url.endswith(".zip"):
            with zipfile.ZipFile(io.BytesIO(data)) as zip_ref:
                member = "mmseqs/bin/mmseqs.exe"
                zip_ref.extract(member, path=install_dir)
                extracted_path = os.path.join(install_dir, member)
                dest_path = os.path.join(install_dir, binary_name)
                shutil.move(extracted_path, dest_path)
    except (KeyError, tarfile.TarError, zipfile.BadZipFile) as e:
        raise RuntimeError(f"Failed to extract MMseqs2 archive: {e}") from e

class CustomBuildCommand(build_py):
    """Custom build command to download MMseqs2 before building."""
    def run(self):
        install_dir = os.path.join(self.build_lib, "pymmseqs", "bin")
        os.makedirs(install_dir, exist_ok=True)
        
        url, binary_name = get_mmseqs_download_info()
        data = download_mmseqs_binary(url)
        extract_mmseqs_binary(data, url, install_dir, binary_name)
        
        super().run()

setup(
    name="pymmseqs",
    version="0.1.0",
    packages=find_packages(),
    package_data={
        "pymmseqs": ["bin/*"]
    },
    include_package_data=True,
    cmdclass={
        "build_py": CustomBuildCommand
    }
)