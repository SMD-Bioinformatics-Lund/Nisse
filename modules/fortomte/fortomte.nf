process CREATE_PEDIGREE_FILE {
    tag "pedigree"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.8.3'
        : 'biocontainers/python:3.8.3'}"

    input:
    val samples

    output:
    path ("*.ped"), emit: ped
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def case_name = samples[0].case
    outfile_text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\\t')
    samples
        .unique { it.sample }
        .each { sample ->
            outfile_text += "\\n" + [sample.case, sample.sample, "0", "0", sample.sex, "affected"].join('\\t')
        }

    // FIXME: Cleanup
    // Original code
    // for(int i = 0; i<samples.size(); i++) {
    //     def sample_name =  samples[i].sample
    //     if (!samples_list.contains(sample_name)) {
    //         outfile_text += "\\n" + [samples[i].case, sample_name, samples[i].paternal, samples[i].maternal, samples[i].sex, samples[i].phenotype].join('\\t')
    //         samples_list.add(sample_name)
    //     }
    // }
    """
    echo -e "${outfile_text}" >${case_name}.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_pedigree_file: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def case_name = samples[0].case
    """
    touch ${case_name}.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_pedigree_file: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


process PEDDY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/peddy:0.4.8--pyh5e36f6f_0'
        : 'biocontainers/peddy:0.4.8--pyh5e36f6f_0'}"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    path ped

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.csv"), emit: csv
    tuple val(meta), path("*.peddy.ped"), emit: ped
    tuple val(meta), path("*.png"), optional: true, emit: png
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    peddy \\
        ${args} \\
        --prefix ${prefix} \\
        --plot \\
        -p ${task.cpus} \\
        ${vcf} \\
        ${ped}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peddy: \$( peddy --version 2>&1 | tail -1 | sed 's/peddy, version //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ped_check.csv
    touch ${prefix}.vs.html
    touch ${prefix}.het_check.csv
    touch ${prefix}.sex_check.csv
    touch ${prefix}.peddy.ped
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peddy: \$( peddy --version 2>&1 | tail -1 | sed 's/peddy, version //' )
    END_VERSIONS
    """
}


process ESTIMATE_HB_PERC {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.base}"

    input:
    tuple val(meta), path(gene_counts)
    path hb_genes

    output:
    path ("${meta.sample}_perc_mapping.json"), emit: tsv

    script:
    """
    calculate_perc_mapping.py \\
        --target_genes "${hb_genes}" \\
        --gene_counts "${gene_counts}" \\
        --strandedness "${meta.strandedness}" \\
        --out_json "${meta.sample}_perc_mapping.json"
    """

    stub:
    """
    touch "${meta.sample}_perc_mapping.json"
    """
}
