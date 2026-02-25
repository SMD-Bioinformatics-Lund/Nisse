process OUTPUT_VERSIONS {

    label "process_low"
    container "${params.containers.base}"

    input:
    val versions_list
    val case_id

    output:
    path ("${case_id}_versions.yaml"), emit: yaml

    script:
    versions_joined = versions_list.join('\n')
    """
    cat <<-END_VERSIONS > "${case_id}_versions.yaml"
    ${versions_joined}
    END_VERSIONS
    """

    stub:
    """
    touch "${case_id}_versions.yaml"
    """
}
