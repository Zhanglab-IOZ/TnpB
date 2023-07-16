This repository holds the pipelines used in the [TnpB project](#).

## Content

The pipelines are organized under subfolders.

- `tam_depletion_seq/`: Identify the depleted TAM from sequencing reads.
- `is605_annotation/`: Annotate *de novo* IS605 elements from prokaryotic genomes.

## Dependencies

- [snakemake](https://snakemake.github.io/) is needed for running the pipelines, and
- [conda](https://conda.io/) or [mamba](https://github.com/mamba-org/mamba) is needed to manage other software dependencies.

The pipelines are developed and tested under snakemake `v7.14.0`.

## Usage

To use these pipelines, first clone this repository and install the dependencies.
Detailed information on preparing input data and running the pipelines is available under the corresponding folder.

## Citation

To cite this repo in publications, please use

> Xiang, G., Li, Y., Sun, J. *et al.* Evolutionary mining and functional characterization of TnpB nucleases identify efficient miniature genome editors. *Nat Biotechnol* (2023). https://doi.org/10.1038/s41587-023-01857-x
