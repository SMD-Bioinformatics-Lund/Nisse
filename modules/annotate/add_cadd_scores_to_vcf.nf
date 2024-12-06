process ADD_CADD_SCORES_TO_VCF {

    tag "${meta.sample}"
    label 'process_low'
    container "${params.containers.genmod}"

    input: 
        tuple val(meta), file(vcf), file(vcf_tbi), file(cadd_scores), file(cadd_scores_tbi)

    output:
        tuple val(meta), file("${meta.sample}.cadd.vcf"), emit: vcf
        tuple val(meta), file("*versions.yml"), emit: versions

    script:
        """
        genmod annotate --cadd-file ${cadd_scores} ${vcf} > ${meta.sample}.cadd.vcf

        ${add_cadd_scores_to_vcf_version(task)}
        """

    stub:
        """
        touch "${meta.sample}.cadd.vcf"
        ${add_cadd_scores_to_vcf_version(task)}
        """
}

def add_cadd_scores_to_vcf_version(task) {
    """
    cat <<-END_VERSIONS > ${task.process}_versions.yml
    ${task.process}:
        genmod: \$(echo \$(genmod --version 2>&1) | sed -e "s/^.*genmod version: //")
    END_VERSIONS
    """
}
