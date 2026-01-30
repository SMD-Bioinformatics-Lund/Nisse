process OUTPUT_VERSIONS {

    label "process_low"
    container "${params.containers.base}"

    input:
    val versions_list

    output:
    path ("versions.yaml"), emit: yaml

    script:
    versions_joined = versions_list.join('\n')
    """
    cat <<-END_VERSIONS > "versions.yaml"
    ${versions_joined}
    END_VERSIONS
    """

    stub:
    """
    touch "versions.yaml"
    """
}
