# ntEdit+Sealer Assembly Finishing Protocol

An automated protocol for finishing long-read genome assemblies using short reads. [ntEdit](https://github.com/bcgsc/ntEdit) polishes the draft assembly and flags erroneous regions, then [Sealer](https://github.com/bcgsc/abyss/tree/master/Sealer) fills assembly gaps and erroneous sequence regions flagged by ntEdit. The protocol is implemented as a Makefile pipeline.

![ntEdit+Sealer protocol flowchart](/ntEdit_Sealer_flowchart.jpg)

## Dependencies

- GNU Make
- Python 3
- [ntHits](https://github.com/bcgsc/nthits) v0.0.1+
- [ntEdit](https://github.com/bcgsc/ntEdit) v1.3.5+
- [ABySS](https://github.com/bcgsc/abyss) v2.3.2+ (includes Sealer and ABySS-Bloom)

## Installation

The ntEdit+Sealer dependencies are available from [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html):

```bash
conda install -c bioconda nthits ntedit abyss
```

All dependencies are also available from [Homebrew](https://docs.brew.sh/Installation):

```bash
brew install brewsci/bio/nthits ntedit abyss
```

This repository, containing the Makefile pipeline and additional scripts, can be cloned from Github:

```bash
git clone https://github.com/bcgsc/ntedit_sealer_protocol.git
```

To run the protocol, ensure that all dependencies are on your `PATH`.

## Example Command

For example, to run the pipeline on a draft long-read assembly `draft-assembly.fa` with short read files `reads_1.fq.gz` and `reads_2.fq.gz`, k-mer lengths `k=80`, `k=65` and `k=50`, specifying the ABySS-Bloom Bloom filter size to be `5G`:

```bash
ntedit-sealer finish seqs=draft-assembly.fa reads='reads_1.fq.gz reads_2.fq.gz' k='80 65 50' b=5G
```

The corrected, finished assembly can be found with the suffix `.ntedit_edited.prepd.sealer_scaffold.fa`.

## Help Page
```
Usage: ntedit-sealer finish [OPTION=VALUE]

General options:
seqs			Draft assembly name [seqs]. File must have .fa extension
reads			Read file(s). All files must have .fq.gz extension. Must be separated by spaces and surrounded by quotes
k			K-mer sizes. List must be descending, separated by spaces and surrounded by quotes
t			Number of threads [8]
time			If True, will log the time for each step [False]

ntEdit options:
X			Ratio of number of kmers in the k subset that should be missing in order to attempt fix (higher=stringent) [0.5]
Y			Ratio of number of kmers in the k subset that should be present to accept an edit (higher=stringent) [0.5]

ABySS-bloom options:
b			Bloom filter size (e.g. 100M)

Sealer options:
L			Length of flanks to be used as pseudoreads [100]
P			Maximum alternate paths to merge; use 'nolimit' for no limit [10]

Notes:
 - Pass all parameter list values (reads, k) as space-separated values surrounded by quotation marks, e.g. k='80 65 50'
 - Ensure that all input files are in the current working directory, making soft-links if needed
 - K-mer lengths will be used in the order they are provided. Ensure that they are sorted in descending order (largest to smallest)
```

Running `ntedit-sealer help` prints the help documentation.
