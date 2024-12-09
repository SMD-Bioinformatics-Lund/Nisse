process OUTPUT_VERSIONS {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.base}"

    input:
    tuple val(meta), path(versions)

    output:
    tuple val(meta), path("versions.yaml"), emit: yaml

    script:
    // versions_joined = versions.join(' ')
    """
    cat ${versions} > "versions.yaml"
    """

    stub:
    """
    touch "versions.yaml"
    """
}
