<div align="center">
<h1>
    PyMMseqs 🚀
</h1>

![GitHub Actions](https://img.shields.io/github/actions/workflow/status/heispv/pymmseqs/pypi-publish.yaml?style=plastic&logo=github-actions&label=CI)
![License](https://img.shields.io/github/license/heispv/pymmseqs?style=plastic&color=orange&logo=github&label=License)
![GitHub stars](https://img.shields.io/github/stars/heispv/pymmseqs?style=social&label=Stars)

PyMMseqs is a powerful Python wrapper for [MMseqs2](https://github.com/soedinglab/MMseqs2). It seamlessly integrates MMseqs2’s advanced functionality into your Python workflows, allowing you to effortlessly execute MMseqs2 commands and parse their outputs into convenient Python objects for further analysis. Whether you're clustering sequences, searching databases, or analyzing large-scale biological data, PyMMseqs simplifies the process while maintaining the performance and flexibility of MMseqs2.
</div>

---

## 🗝️ Features

- **Seamless Integration**: Execute MMseqs2 commands directly within your Python code, eliminating the need for shell scripting or external command-line tools.
- **Output Parsing**: Convert MMseqs2 outputs into Python objects (e.g., generators, dictionaries) for easy manipulation and analysis.
- **High Performance**: Leverage the speed and efficiency of MMseqs2 while enjoying the flexibility of Python.
- **Cross-Platform**: Use PyMMseqs via pip or Docker, ensuring compatibility across different environments.

---

## 🛠️ Installation

PyMMseqs can be installed in two ways: via pip (recommended for most users) or using a Docker image (ideal for reproducible environments).

### Installing via pip
The `pymmseqs` package is currently available on TestPyPI. To install it, use the following command:

```bash
pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple pymmseqs
```

> **Note**: This command uses TestPyPI as the primary source for `pymmseqs` while fetching dependencies from the main PyPI index. Once the package is available on PyPI, this command will be simplified to `pip install pymmseqs`.

> **Important**: All dependencies, including MMseqs2, are automatically installed and configured when you install `pymmseqs` via pip. No additional setup is required.

### Using Docker Image
For users who prefer not to install PyMMseqs locally or want a pre-configured environment, a Docker image is available on GitHub Container Registry (GHCR).

#### Debian-based Image
To pull the Debian-based Docker image, run:

```bash
docker pull ghcr.io/heispv/pymmseqs:debian
```

> **Tip**: Using Docker ensures that all dependencies, including MMseqs2, are pre-installed and configured, making it ideal for reproducible workflows.

---

## 🚀 Quick Start

Here’s a simple example to get you started with PyMMseqs. This example demonstrates how to perform sequence clustering and parse the results.

If you were using MMseqs2 directly in the terminal, you would run the following command to cluster sequences:

```bash
mmseqs easy-cluster human.fasta human_clust tmp --min-seq-id 0.9
```

With PyMMseqs, you can achieve the same result directly in Python, and parse the output to Python objects for further analysis.

```python
from pymmseqs.commands import easy_cluster

# Perform clustering on a FASTA file (equivalent to the terminal command above)
human_cluster = easy_cluster("human.fasta", "human_clust", "tmp", min_seq_id=0.9)

# Convert the clustering output to a Python generator for efficient iteration
cluster_gen = human_cluster.to_gen()

# Analyze the results: Find and print the representative sequence of a large cluster (>100 members)
for cluster in cluster_gen:
    if len(cluster["members"]) > 100:
        print(f"Representative sequence of a large cluster: {cluster['rep']}")
        break
```
---

## 📖 Documentation

For detailed usage instructions, advanced examples, and API references, please visit the [PyMMseqs Wiki](https://github.com/heispv/pymmseqs/wiki).

---

## 🔧 Prerequisites

To use PyMMseqs, you only need:
- **Python**: Version 3.10 or higher.

> **Note**: All other dependencies, including MMseqs2, are automatically installed when you install `pymmseqs` via pip or use the Docker image.

---

## 🤝 Contributing

We welcome contributions to PyMMseqs! If you’d like to contribute, please:
1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/YourFeature`).
3. Commit your changes (`git commit -m "Add YourFeature"`).
4. Push to the branch (`git push origin feature/YourFeature`).
5. Open a Pull Request.

For bug reports, feature requests, or questions, please open an issue on the [GitHub Issues page](https://github.com/heispv/pymmseqs/issues).

---

## 📜 License

PyMMseqs is licensed under the [MIT License](LICENSE). See the [LICENSE](LICENSE) file for more details.

---

## 🌟 Support

If you find PyMMseqs useful, please consider giving the repository a star on GitHub! ⭐ It helps others discover the project and motivates further development.

For questions, feedback, or support, feel free to open an issue or contact the maintainers.
