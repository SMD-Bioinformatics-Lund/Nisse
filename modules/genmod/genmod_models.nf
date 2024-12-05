process GENMOD_MODELS {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.genmod}"

    input:
        tuple val(meta), path(vcf)
        tuple val(meta2), path(family)

    output:
        tuple val(meta), path("${meta.sample}_models.vcf"), emit: vcf
        path("versions.yaml"), emit: versions

    script:
    """
    genmod \\
        models \\
        --whole_gene \\
        --processes "${task.cpus}" \\
        --family_file "${family}" \\
        --outfile "${meta.sample}_models.vcf" \\
        "${vcf}"
    
    echo "FIXME" > "versions.yaml"
    """

    stub:
    """
    touch "${meta.sample}_models.vcf"
    echo "FIXME" > "versions.yaml"
    """
}