# setup.py

import os
import platform
import tarfile
import zipfile
import io
from urllib.request import urlopen
from urllib.error import URLError
from setuptools import setup, find_packages
import shutil
import sysconfig
import glob


def download_mmseqs2():
    """Download and extract the mmseqs2 binary."""
    print("Starting mmseqs2 download process...")

    # get_python_lib = sysconfig.get_path

    # Define the directory where mmseqs2 will be stored
    install_dir = os.path.join(os.path.dirname(__file__), 'pymseqs', 'bin')
    # temp_dir = os.path.join(get_python_lib('purelib'), 'pymseqs', 'temp')
    os.makedirs(install_dir, exist_ok=True)
    # os.makedirs(temp_dir, exist_ok=True)
    print(f"Installation directory for mmseqs2: {install_dir}")

    # Determine the OS and architecture
    system = platform.system()
    machine = platform.machine().lower()
    print(f"Detected system: {system}, architecture: {machine}")

    # Define binary name and URL based on OS and architecture
    if system == 'Linux':
        if 'x86_64' in machine or 'amd64' in machine:
            binary_name = 'mmseqs'
            url = 'https://github.com/soedinglab/mmseqs2/releases/download/v13.5/mmseqs-linux-amd64.tar.gz'
        elif 'arm' in machine or 'aarch64' in machine:
            binary_name = 'mmseqs'
            url = 'https://github.com/soedinglab/mmseqs2/releases/download/v13.5/mmseqs-linux-arm64.tar.gz'
        else:
            raise RuntimeError(f"Unsupported architecture: {machine}")
    elif system == 'Darwin':
        if 'x86_64' in machine:
            binary_name = 'mmseqs'
            url = 'https://github.com/soedinglab/mmseqs2/releases/download/v13.5/mmseqs-osx-amd64.tar.gz'
        elif 'arm64' in machine or 'aarch64' in machine:
            binary_name = 'mmseqs'
            url = 'https://github.com/soedinglab/MMseqs2/releases/download/16-747c6/mmseqs-osx-universal.tar.gz'
        else:
            raise RuntimeError(f"Unsupported architecture: {machine}")
    elif system == 'Windows':
        if 'x86_64' in machine or 'amd64' in machine:
            binary_name = 'mmseqs.exe'
            url = 'https://github.com/soedinglab/mmseqs2/releases/download/v13.5/mmseqs-win64.zip'
        else:
            raise RuntimeError(f"Unsupported architecture: {machine}")
    else:
        raise RuntimeError(f"Unsupported operating system: {system}")

    # Define the path to the binary
    binary_path = os.path.join(install_dir, binary_name)

    # Check if the binary already exists to avoid re-downloading
    if os.path.exists(binary_path):
        print(f"mmseqs2 binary already exists at {binary_path}")
        return

    print(f"Downloading mmseqs2 from {url}...")

    try:
        response = urlopen(url)
        data = response.read()
        print(f"Downloaded data size: {len(data)} bytes")
    except URLError as e:
        raise RuntimeError(f"Failed to download mmseqs2 from {url}: {e}")

    # Extract the binary based on the file type
    if url.endswith('.tar.gz'):
        try:
            with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
                print(f"Extracting mmseqs from tar.gz archive...")
                # Extract the binary directly to the desired directory
                member = tar.getmember('mmseqs/bin/mmseqs') if system != 'Windows' else tar.getmember('mmseqs/bin/mmseqs.exe')
                member.name = os.path.basename(member.name)  # Flatten the path
                tar.extract(member, path=install_dir)
                binary_path = os.path.join(install_dir, binary_name)
                print(f"Binary extracted to {binary_path}")
        except KeyError:
            raise RuntimeError(f"Binary 'mmseqs/bin/{binary_name}' not found in the tar.gz archive.")
        except Exception as e:
            raise RuntimeError(f"Failed to extract mmseqs2 from tar.gz: {e}")
    elif url.endswith('.zip'):
        try:
            with zipfile.ZipFile(io.BytesIO(data)) as zip_ref:
                print(f"Extracting mmseqs.exe from zip archive...")
                # Extract the binary directly to the desired directory
                member = 'mmseqs/bin/mmseqs.exe'
                zip_ref.extract(member, path=install_dir)
                extracted_path = os.path.join(install_dir, 'mmseqs', 'bin', binary_name)
                shutil.move(extracted_path, binary_path)  # Adjust path to move to the correct directory
                print(f"Moved binary from {extracted_path} to {binary_path}")
        except KeyError:
            raise RuntimeError(f"Binary 'mmseqs/bin/{binary_name}' not found in the zip archive.")
        except Exception as e:
            raise RuntimeError(f"Failed to extract mmseqs2 from zip: {e}")

        



# Run the binary download process before setup
download_mmseqs2()

setup(
    name='pymseqs',
    version='0.1.0',
    description='Python wrapper for mmseqs2',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='heispv',
    author_email='peymanvahidi1998@gmail.com',
    url='https://github.com/heispv/pymseqs',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'pymseqs': ['bin/*'],
    },
    # install_requires=[
    #     Add your package dependencies here
    # ],
    entry_points={
        'console_scripts': [
            'pymseqs=pymseqs:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)

artifacts = ['build', 'dist', '*.egg-info', 'pymseqs/bin']
current_dir = os.getcwd()  # Get the current working directory
system = platform.system()


for pattern in artifacts:
    # Prepend current_dir to each artifact and use glob to find matches
    for path in glob.glob(os.path.join(current_dir, pattern)):
        try:
            if os.path.isdir(path):
                shutil.rmtree(path)
            else:
                os.remove(path)
            print(f"Removed {path}")
        except Exception as e:
            print(f"Failed to remove {path}: {e}")