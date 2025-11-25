Home for CMD [Tomte](https://github.com/genomic-medicine-sweden/tomte) post processing pipeline.

### Setting things up

(With GitHub access, instructions for hopper below). Retrieving the correct versions of Nisse, Tomte and config-files.

```
$ git clone --recurse-submodules git@github.com:SMD-Bioinformatics-Lund/Nisse.git
```

Note that this will not work on hopper due to lack of internet access. There, you will have to do something like the following:

```
$ git clone /fs1/pipelines/bare/nisse
$ cd nisse
$ git config -f .gitmodules submodule.config-files.url /fs1/pipelines/bare/config-files
$ git config -f .gitmodules submodule.tomte.url        /fs1/pipelines/bare/tomte
# Confirm the updated paths
$ cat .gitmodules
[submodule "tomte"]
        path = tomte
        url = /fs1/pipelines/bare/tomte
[submodule "config-files"]
        path = config-files
        url = /fs1/pipelines/bare/config-files
# FIXME: What is happening in these steps
$ git submodule sync --recursive
$ git submodule update --init --recursive
```

Furthermore, certain Tomte files will need to be linked into the Nisse repo.

```
cd Nisse
ln -s tomte/nextflow_schema.json
ln -s tomte/assets
```

And scripts from both Nisse (`bin_nisse`) and Tomte (`tomte/bin`) needs to be copied into `bin`.

```
mkdir bin/
cp bin_nisse/* bin/
cp tomte/bin/* bin/
```

### Running Nisse

Run Tomte as part of Nisse. Tomte and config-files repo are available as sub modules.

Execute while providing the Tomte configs:

```
nextflow run main.nf \
    --csv ... \
    --outdir ... \
    -c ./tomte/nextflow_tomte.config \
    -c ./config-files/nextflow/cmd_tomte.config \
    -c ./config-files/nextflow/cmd_nisse.config
```

### Running Nisse

Takes the CSV-samplesheet used to run Tomte and and path to the Tomte results folder as input.

An example input CSV using HG002 GIAB downsampled to 15M can be found at `/fs2/resources/ref/hg38/tomte/testdata/hg002.csv`.

Example run script which can be executed using `sbatch`:

```
#!/bin/bash
#SBATCH --job-name=nisse_test_run
#SBATCH --output=%j.log
#SBATCH --ntasks=4
#SBATCH --mem=4gb
#SBATCH --time=7-00:00:00
#SBATCH --partition="grace-lowest"

module load Java/13.0.2
module load nextflow/24.04.3
module load singularity/3.2.0

NXF_OFFLINE=true

timestamp=$(date +"%y%m%d_%H%M")
outdir="/nisse/results/${timestamp}"
mkdir -p ${outdir}

tomte_base_config="./tomte/nextflow_tomte.config"
tomte_smd_config="./config-files/nextflow/cmd_tomte.config"
nisse_smd_config="./config-files/nextflow/cmd_nisse.config"

csv="/path/to/samplesheet.csv"
work="/path/to/workdir"
main_nf="/base/nisse/main.nf"

nextflow run "${main_nf}" \
    -c "${tomte_base_config}" \
    -c "${tomte_smd_config}" \
    -c "${nisse_smd_config}" \
    --outdir "${outdir}" \
    -work-dir "${work}" \
    --input "${csv}" | tee "${outdir}/run.log"
```

### Map HGNC symbol to ID

A two-column file mapping HGNC symbol to ID is currently needed, as this isn't performed by Tomte but required by Scout (`params.hgnc_map`).

Steps to prepare:

1. Download the latest tab `hgnc_complete_set` file from [https://www.genenames.org/download/archive/#!/#tocAnchor-1-3](genenames.org)
2. Extract the symbol and ID columns, without header

```
tail -n +2 hgnc_complete_set.txt | head | cut -f1,2 | sed "s/HGNC://" > hgnc_id_symbol_map.tsv
```

