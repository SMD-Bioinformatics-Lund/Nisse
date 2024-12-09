process OUTPUT_VERSIONS {

    label "process_low"
    container "${params.containers.base}"

    input:
    path versions

    output:
    path ("versions.yaml"), emit: yaml

    script:
    versions_joined = versions.join(' ')
    """
    cat ${versions_joined} > "versions.yaml"
    """

    stub:
    """
    touch "versions.yaml"
    """
}
