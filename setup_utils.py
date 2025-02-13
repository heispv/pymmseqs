# setup_utils.py
import io
import os
import platform
import shutil
import tarfile
import zipfile
from urllib.error import URLError
from urllib.request import urlopen

MMSEQS_VERSION = "16-747c6"

def get_mmseqs_download_info():
    """Determine the appropriate MMseqs2 binary URL and filename for the current platform."""
    system = platform.system()
    machine = platform.machine().lower()

    base_url = "https://github.com/soedinglab/mmseqs2/releases/download"

    if system == "Linux":
        if machine in ("x86_64", "amd64"):
            return f"{base_url}/{MMSEQS_VERSION}/mmseqs-linux-avx2.tar.gz", "mmseqs"
        if machine in ("arm", "aarch64"):
            return f"{base_url}/{MMSEQS_VERSION}/mmseqs-linux-arm64.tar.gz", "mmseqs"
    elif system == "Darwin":
        if machine in ("x86_64", "amd64"):
            return f"{base_url}/{MMSEQS_VERSION}/mmseqs-osx-arm64.tar.gz", "mmseqs"
        if machine in ("arm64", "aarch64"):
            return f"{base_url}/{MMSEQS_VERSION}/mmseqs-osx-universal.tar.gz", "mmseqs"
    elif system == "Windows":
        if machine in ("x86_64", "amd64"):
            return f"{base_url}/{MMSEQS_VERSION}/mmseqs-win64.zip", "mmseqs.exe"

    raise RuntimeError(f"Unsupported platform: {system} {machine}")

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
                member_path = "mmseqs/bin/mmseqs.exe" if platform.system() == "Windows" else "mmseqs/bin/mmseqs"
                member = tar.getmember(member_path)
                member.name = os.path.basename(member.name)  # Flatten directory structure
                tar.extract(member, path=install_dir)
                print(f"Extracted {member.name} to {install_dir}")
        elif url.endswith(".zip"):
            with zipfile.ZipFile(io.BytesIO(data)) as zip_ref:
                member = "mmseqs/bin/mmseqs.exe"
                zip_ref.extract(member, path=install_dir)
                extracted_path = os.path.join(install_dir, member)
                dest_path = os.path.join(install_dir, binary_name)
                shutil.move(extracted_path, dest_path)
                print(f"Moved {member} to {dest_path}")
    except (KeyError, tarfile.TarError, zipfile.BadZipFile) as e:
        raise RuntimeError(f"Failed to extract MMseqs2 archive: {e}") from e

def download_mmseqs2(install_dir):
    """Download and install MMseqs2 binary to the specified directory."""
    print("Starting MMseqs2 installation process...")
    os.makedirs(install_dir, exist_ok=True)
    print(f"Installation directory: {install_dir}")

    binary_name = "mmseqs.exe" if platform.system() == "Windows" else "mmseqs"
    binary_path = os.path.join(install_dir, binary_name)
    
    if os.path.exists(binary_path):
        print(f"MMseqs2 already exists at {binary_path}")
        return

    url, _ = get_mmseqs_download_info()
    data = download_mmseqs_binary(url)
    extract_mmseqs_binary(data, url, install_dir, binary_name)
