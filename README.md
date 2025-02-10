Home for CMD [Tomte](https://github.com/genomic-medicine-sweden/tomte) post processing pipeline.

### Running Nisse (new)

Run Tomte as part of Nisse. Requires that you have:

* The Tomte repository cloned into the Nisse as a subfolder named "tomte"

Execute while providing the Tomte configs:

```
nextflow run main.nf \
    --input ... \
    --outdir ... \
    -c nextflow.config \
    -c /fs1/jakob/src/config-files/nextflow/cmd_tomte.config \
    -c tomte/nextflow.config
```

### Running Nisse

Takes the CSV-samplesheet used to run Tomte and and path to the Tomte results folder as input.

Example run script which can be executed using `sbatch`:

```
#!/bin/bash
#SBATCH --job-name=nisse
#SBATCH --output=slurm_logs/%j.log
#SBATCH --ntasks=4
#SBATCH --mem=4gb
#SBATCH --time=7-00:00:00
#SBATCH --partition="grace-lowest"

module load Java/13.0.2
module load nextflow/24.04.3
module load singularity/3.2.0

export NXF_HOME="/local/nextflowconf/jakob/"
export NXF_PLUGINS_DIR="/local/nextflowconf/jakob/plugins/"

timestamp=$(date +"%y%m%d_%H%M")
outdir="/nisse/results/${timestamp}"
mkdir -p ${outdir}

csv="/path/to/samplesheet.csv"
tomte_results="/path/to/tomte_results"
work="/path/to/workdir"
main_nf="/base/nisse/main.nf"

nextflow run "${main_nf}" \
    -work-dir "${work}" \
    --tomte_results "${tomte_results}" \
    --input "${csv}" \
    --outdir ${outdir} | tee "${outdir}/run.log"
```

### Map HGNC symbol to ID

A two-column file mapping HGNC symbol to ID is currently needed, as this isn't performed by Tomte but required by Scout (`params.hgnc_map`).

Steps to prepare:

1. Download the latest tab `hgnc_complete_set` file from [https://www.genenames.org/download/archive/#!/#tocAnchor-1-3](genenames.org)
2. Extract the symbol and ID columns, without header

```
tail -n +2 hgnc_complete_set.txt | head | cut -f1,2 | sed "s/HGNC://" > hgnc_id_symbol_map.tsv
```

