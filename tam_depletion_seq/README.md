# tam_depletion_seq

```mermaid
graph LR
    0(["all"]);
    1(["plot_logo"]);
    2(["build_logodata"]);
    3(["test_depletion"]);
    4(["merge_pairs"]);
    5(["find_tam"]);
    1 --> 0
    2 --> 1
    3 --> 2
    4 --> 3
    5 --> 4
```

Identify depleted TAM from seqeuencing data and build the sequence logo.

## Quick start guide

### Get the workflow

Clone this repo and set the working directory:

```bash
git clone https://github.com/Zhanglab-IOZ/TnpB.git
cd TnpB/tam_depletion_seq
```

### Prepare sequencing data

```bash
mkdir -p resources
```

Download the following files to the `resources/` folder:

- `TAM_NC_R1.fq.gz`, `TAM_NC_R2.fq.gz`: The control group.
- `TAM_ISHahl1_R1.fq.gz`, `TAM_ISHahl1_R2.fq.gz` and `TAM_ISTfu1_R1.fq.gz`, `TAM_ISTfu1_R2.fq.gz`: The 2 test groups.

### Create config file

Create the config file config/config.yaml using the template file:

```bash
cp config/config_example.yaml config/config.yaml
```

Edit the config file accordingly if different data files are used.

### Execute the workflow

Run the workflow with snakemake:

```bash
snakemake --use-conda --cores all
```

### Result files

The resulting TAM logo is avaliable as `results/logo/*.logo.svg`.
