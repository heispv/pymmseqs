# PyMMseqs

A Python wrapper for [MMseqs2](https://github.com/soedinglab/MMseqs2), designed to streamline the integration of MMseqs2 functionality into your Python workflows. 🚀

## Installation

### Installing via pip
The `pymmseqs` package is currently available on TestPyPI. You can install it using pip with the following command:

```bash
pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple pymmseqs
```

This command installs `pymmseqs` from TestPyPI while also using the main PyPI index to fetch any necessary dependencies.

## Simple Example

```python
from pymmseqs.commands import easy_cluster

# Perform clustering
cluster = easy_cluster(
    fasta_files = "human.fasta"
    cluster_prefix = "output/human_clust",
    tmp_dir = "output/tmp",
    min_seq_id = 0.9,
    c: float = 0.7
)

# Parse output to python generator
cluster_gen = cluster.to_gen()

# Print a cluster representative of a cluster which has more than 100 members
for cluster in cluster_gen:
    if len(cluster["members"]) > 100:
        print(cluster["rep"])
        break

```

## Documentation
More detailed documentation can be found in [Wiki](https://github.com/heispv/pymmseqs/wiki)
