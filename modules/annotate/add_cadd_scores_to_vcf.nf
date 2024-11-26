process ADD_CADD_SCORES_TO_VCF {
    cpus 4
    tag "$group"
    memory '1 GB'
    time '5m'
    container "${params.container_genmod}"

    input: 
        set group, file(vcf), file(cadd_scores), file(cadd_scores_tbi)

    output:
        set group, file("${group}.cadd.vcf")
        set group, file("*versions.yml")

    script:
        """
        genmod annotate --cadd-file ${cadd_scores} ${vcf} > ${group}.cadd.vcf

        ${add_cadd_scores_to_vcf_version(task)}
        """

    stub:
        """
        touch "${group}.cadd.vcf"
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
