process GENMOD_COMPOUND {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.genmod}"

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("*_compound.vcf"), emit: vcf
        path("versions.yaml"), emit: versions

    script:
    """
    genmod \\
        compound \\
        --processes ${task.cpus} \\
        --outfile ${meta.sample}_compound.vcf \\
        --penalty 0 \\
        ${vcf}
    echo "FIXME" > "versions.yaml"
    """

    stub:
    """
    touch ${meta.sample}_compound.vcf
    echo "FIXME" > "versions.yaml"
    """
}