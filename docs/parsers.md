# pymmseqs.parsers
This module is for parsing the output of MMseqs2 commands as Python objects, making it easier to work with the results in your Python scripts.

## Important Things to Know

- You obtain a parser object when using the `pymmseqs.commands` module to execute a command.
- Pro users can pass the `pymmseqs.config` object to the parser after execution to leverage additional features of `pymmseqs.parsers`.
- Refer to the sections below for detailed instructions on how to use parser objects effectively.

---

## [CreateDBParser](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/parsers/createdb_parser.py)
This parser processes the output of the MMseqs2 `createdb` command.

### For Normal Users
When using the `createdb` from `pymmseqs.commands`, you receive a `CreateDBParser` object.

Example:
```python
from pymmseqs.commands import createdb

query_database = createdb(
    fasta_file="data/query.fasta",
    sequence_db="output/query_db",
)

# Retrieve the path to the database
print(query_database.get_path())
```

### For Pro Users
Pro users can utilize the `pymmseqs.config.CreateDBConfig` object for additional flexibility.

Example:
```python
from pymmseqs.config import CreateDBConfig
from pymmseqs.parsers import CreateDBParser

# Create the configuration
createdb_config = CreateDBConfig(
    fasta_file="data/query.fasta",
    sequence_db="output/query_db",
    dbtype=1,
)

# Execute the configuration
createdb_config.run()

# Obtain the parser object
createdb_parser = CreateDBParser(createdb_config)

# Retrieve the path to the database
print(createdb_parser.get_path())
```

**Want to know what happens when you run `config.run()`?**

Check out the [`pymmseqs.config`](../pymmseqs.config) module for details.

---

## [EasyClusterParser](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/parsers/easy_cluster_parser.py)
The `EasyClusterParser` provides several methods for handling clustering results:

### Methods:

#### `to_list()`
- Returns a list of dictionaries, each representing a cluster.
- Each dictionary contains:
    - **`rep`**: The representative sequence ID.
    - **`members`**: A list of dictionaries for cluster members with:
        - `seq_id`: Unique sequence identifier.
        - `header`: Full FASTA header.
        - `sequence`: Sequence data (nucleotide or protein).

**Example Output:**
```python
[
    {
        "rep": "seq1",
        "members": [
            {"seq_id": "seq1", "header": "sp|seq1", "sequence": "SEQWENCE"},
            {"seq_id": "seq2", "header": "sp|seq2", "sequence": "PRTEINSEQWENCE"},
        ]
    },
    {
        "rep": "seq3",
        "members": [
            {"seq_id": "seq3", "header": "sp|seq3", "sequence": "PRTEIN"}
        ]
    }
]
```

#### `to_pandas()`
- Converts cluster data into a **pandas DataFrame** for easier analysis.
- Columns include:
    - `rep`: Representative sequence ID (index)
    - `seq_id`: Sequence identifier
    - `header`: FASTA header
    - `sequence`: Sequence data

#### `to_gen()`
- Returns a **generator** of clusters, allowing memory-efficient iteration.

#### `to_path()`
- Returns a list of relevant output file paths:
    - `cluster_prefix_all_seqs.fasta`: All sequences in clusters
    - `cluster_prefix_cluster.tsv`: Cluster information
    - `cluster_prefix_rep_seqs.fasta`: Representative sequences

**Note:** `to_list()`, `to_pandas()`, and `to_gen()` rely on `cluster_prefix_all_seqs.fasta`.

### For Normal Users
Using the `pymmseqs.commands.easy_cluster` command provides an `EasyClusterParser` object.

Example:
```python
from pymmseqs.commands import easy_cluster

my_cluster = easy_cluster(
    fasta_files="query.fasta",
    cluster_prefix="output/cluster",
    tmp_dir="output/tmp",
    min_seq_id=0.3,
)

# If you have a small dataset...
cluster_list = my_cluster.get_list()
print(cluster_list[:2])
```

### For Pro Users
Pro users can utilize the `pymmseqs.config.EasyClusterConfig` object for additional control.

Example:
```python
from pymmseqs.config import EasyClusterConfig
from pymmseqs.parsers import EasyClusterParser

# Create the configuration
easy_cluster_config = EasyClusterConfig(
    fasta_files="data/query.fasta",
    cluster_prefix="output/cluster",
    tmp_dir="output/tmp",
    min_seq_id=0.3,
)

# Execute the configuration
easy_cluster_config.run()

# Obtain the parser object
easy_cluster_parser = EasyClusterParser(easy_cluster_config)

# Create a generator for clusters
cluster_gen = easy_cluster_parser.get_gen()

# Retrieve the representative sequence of a cluster with more than 10 members
for cluster in cluster_gen:
    if len(cluster["members"]) > 10:
        print(cluster["rep"])
        break
```

---

## [SearchParser](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/parsers/search_parser.py)
The `SearchParser` processes the output of the MMseqs2 `search` command, which performs sequence similarity searches between a query database and a target database.

### Methods:

#### `to_list()`
- Returns a list of dictionaries, each representing a row in the alignment file.
- Each dictionary contains keys corresponding to the alignment data, such as:
  - `query`
  - `target`
  - `fident`
  - `alnlen`
  - `mismatch`
  - `gapopen`
  - `qstart`
  - `qend`
  - `tstart`
  - `tend`
  - `evalue`
  - `bits`

#### `to_pandas()`
- Returns a pandas DataFrame containing the alignment data.
- The columns are the same as the keys in the dictionaries returned by `to_list()`.

#### `to_gen()`
- Returns a generator of dictionaries, each representing a row in the alignment file.
- The keys are the same as the ones in the dictionaries returned by `to_list()`.

#### `to_path()`
- Returns the path to the alignment database.

**Note:** When using `to_list()`, `to_pandas()`, or `to_gen()` methods, the `SearchParser` automatically runs the `mmseqs convertalis` command to convert the binary alignment database to a readable TSV format if needed.

### For Normal Users
When using the `pymmseqs.commands.search` command, you receive a `SearchParser` object.

Example:
```python
from pymmseqs.commands import search
search_result = search(
query_db="output/query_db",
target_db="output/target_db",
    alignment_db="output/search_results",
    tmp_dir="output/tmp",
    min_seq_id=0.3,
    max_seqs=1000,
)

# Get all alignments with sequence identity > 50%
filtered_alignments = []
for alignment in search_result.to_gen():
    if float(alignment["fident"]) > 0.5:
        filtered_alignments.append(alignment)
print(f"Found {len(filtered_alignments)} alignments with >50% identity")
```

### For Pro Users
Pro users can utilize the `pymmseqs.config.SearchConfig` object for additional control.

Example:
```python
from pymmseqs.config import SearchConfig
from pymmseqs.parsers import SearchParser

# Create the configuration with advanced parameters
search_config = SearchConfig(
    query_db="output/query_db",
    target_db="output/target_db",
    alignment_db="output/search_results",
    tmp_dir="output/tmp",
    e=1e-5,
    min_seq_id=0.3,
    c=0.8,
    a=True,
)

# Turn of the logging file saving
search_config.set_config_options(has_log=False)

# Execute the configuration
search_config.run()

# Obtain the parser object
search_parser = SearchParser(search_config)

# Convert results to pandas DataFrame for analysis
results_df = search_parser.to_pandas()

# Filter and analyze results
high_quality_hits = results_df[(results_df['fident'] > 0.7) & (results_df['evalue'] < 1e-10)]
print(f"Found {len(high_quality_hits)} high-quality alignments")
```

---

## [EasySearchParser](https://github.com/heispv/pymmseqs/blob/master/pymmseqs/parsers/easy_search_parser.py)
One of the main differences between the `EasySearchParser` and the `SearchParser` is that the `EasySearchParser` also accepts the inputs as a fasta files, but in the case of the `SearchParser` you need to pass the database paths.

* When we are running the `easy_search` command it will run a `EasySearchConfig` under the hood, with the `format_mode` set to 4. So, we get a .tsv file as output with headers which can be later parsed by the `EasySearchParser`.

### Methods:

#### `to_list()`
- Returns a list of dictionaries, each representing a row in the alignment file.
- Each dictionary contains the following keys (these are the defined in the `format_output` parameter of the `EasySearchConfig` object):
    - `query`
    - `target`
    - `fident`
    - `alnlen`
    - `mismatch`
    - `gapopen`
    - `qstart`
    - `qend`
    - `tstart`
    - `tend`
    - `evalue`
    - `bits`


#### `to_pandas()`
- Returns a pandas DataFrame containing the alignment data.
- The columns are the same as the keys in the dictionaries returned by `to_list()`.

#### `to_gen()`
- Returns a generator of dictionaries, each representing a row in the alignment file.
- The keys are the same as the ones in the dictionaries returned by `to_list()`.

#### `to_path()`
- Returns a list of file paths for the output files.

### Normal Users
When using the `pymmseqs.commands.easy_search` command, you get a `EasySearchParser` object.
```python
from pymmseqs.commands import easy_search

my_search = easy_search(
    query_fasta="query.fasta",
    target_fasta_or_db="target.fasta",
    alignment_file="output/results",
    tmp_dir="output/tmp",
    min_seq_id=0.5,
    max_seqs=10000,
)

```

Since you get a parser object...

Let's say you want to get all the alignments with a sequence identity greater than 50%.
```python
filtered_alignments = []
for alignment in my_search.to_gen():
    if alignment["fident"] > 0.5:
        filtered_alignments.append(alignment)
```

Maybe your dataset is too large to fit into memory...

```python
def get_filtered_alignments(parser, threshold=0.5):
    for alignment in parser.to_gen():
        if alignment["fident"] > threshold:
            yield alignment
```

Then just use the generator to get the filtered alignments, and you don't need to worry about memory issues.

### Pro Users
Pro users can utilize the `pymmseqs.config.EasySearchConfig` object for additional control.

Example:
```python
from pymmseqs.config import EasySearchConfig
from pymmseqs.parsers import EasySearchParser

# Create the configuration
easy_search_config = EasySearchConfig(
    query_fasta="query.fasta",
    target_fasta_or_db="target.fasta",
    alignment_file="output/results",
    tmp_dir="output/tmp",
    realign=True,
    format_mode=4,
)

# Note that if you don't set the format_mode to 4, you will get an error.

# Execute the configuration
easy_search_config.run()

# Obtain the parser object
easy_search_parser = EasySearchParser(easy_search_config)

# Get the alignment list
alignment_list = easy_search_parser.to_list()

```
