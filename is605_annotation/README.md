# is605_annotation workflow

This workflow *de novo* annotate multicopy IS605 elements from prokaryotic genomes.

## Quick start guide

### Get the workflow

To run this workflow, first clone this repo and set the working directory:

```bash
git clone https://github.com/Zhanglab-IOZ/TnpB.git
cd TnpB/is605_annotation
```

### Prepare genome files and the HMM profiles

Download genome assemblies and annotations from NCBI, you may change the accession list as needed:

```bash
mkdir -p resources
cat <<EOF |
GCF_001634285.1_ASM163428v1
GCF_002243625.1_ASM224362v1
GCF_002368015.1_ASM236801v1
GCF_002368335.1_ASM236833v1
GCF_003096215.1_ASM309621v1
GCF_003614235.1_ASM361423v1
GCF_015351765.1_ASM1535176v1
GCF_017829975.1_ASM1782997v1
GCF_017873395.1_ASM1787339v1
GCF_018343715.1_ASM1834371v1
GCF_021390075.1_ASM2139007v1
GCF_900162105.1_PRJEB19297
GCF_902143525.1_P7770
EOF
sed --regexp-extended 's%^(GC.)_([[:digit:]]{3})([[:digit:]]{3})([[:digit:]]{3})%\1/\2/\3/\4\t\0%' |
while read path accession; do
    mkdir -p resources/genomes/${accession}/
    for suffix in genomic.fna genomic.gff translated_cds.faa; do
        wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/${path}/${accession}/${accession}_${suffix}.gz" -O resources/genomes/${accession}/${accession}_${suffix}.gz
    done
done
```

Download the HMM models for TnpB and TnpA from [Pfam](https://www.ebi.ac.uk/interpro/) and [Altae-Tran *et al.* (2021)](https://doi.org/10.1126/science.abj6856), respectively:

```bash
mkdir -p resources/hmm; pushd resources/hmm
wget https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF12323?annotation=hmm -O PF12323.hmm
wget https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF01385?annotation=hmm -O PF01385.hmm
wget https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF07282?annotation=hmm -O PF07282.hmm

wget 'https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.abj6856&file=science.abj6856_Data_S1_to_S4.zip' -O science.abj6856_Data_S1_to_S4.zip
unzip science.abj6856_Data_S1_to_S4.zip
unzip science.abj6856_Data-S1.zip

zcat PF12323.hmm PF01385.hmm PF07282.hmm > ../tnp_models.hmm
cat hhblits_full_TnpA.HMM >> ../tnp_models.hmm
popd; rm -r resources/hmm
```

### Create the config

(Optional) List known TnpB protein sequences in fasta format to filter them out in the result:
```bash
# touch resources/known_tnpb.faa
```

Create the config file `config/config.yaml` using the template file:
```bash
cp config/config_example.yaml config/config.yaml
```
If you use different genome files or want to filter known elements, change the config file accordingly.

### Workflow execution

The workflow contains two parts. First, multicopy IS605 elements are annotated in seperate genomes by `do_genomes`:

```bash
snakemake --use-conda --cores all do_genomes
```

And these elements are clustered to remove redundancy. Also, different copys are aligned to facilitate annotation of element margin. These are performed by `do_elements`:

```bash
snakemake --use-conda --cores all --cores 20 do_elements
```

### Result files

- `results/dedup/candidates.faa` contains a nonredundant set of *de novo* TnpB proteins.

- `results/elements/*/mafft.fna` contains the alignment of each multicopy elements.
Manunal curation is needed for precise margin of the IS element.
