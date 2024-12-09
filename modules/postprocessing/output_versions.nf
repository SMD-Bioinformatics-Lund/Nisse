process OUTPUT_VERSIONS {

    label "process_low"
    container "${params.containers.base}"

    input:
    val versions_list

    output:
    path ("versions.yaml"), emit: yaml

    script:
    versions_joined = versions_list.join(' ')
    """
    cat ${versions_joined} > "versions.yaml"
    """

    stub:
    """
    touch "versions.yaml"
    """
}
