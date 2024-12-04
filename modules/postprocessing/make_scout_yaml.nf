process MAKE_SCOUT_YAML {

    tag "${meta.sample}"

    input:
    tuple val(meta), path(input_file)

    output:
    path 'hello_out.txt', emit: yaml

    script:
    """
    bash hello.sh ${input_file} > "${meta.sample}_out.txt"
    """

    stub:
    """
    touch "${meta.sample}_out.txt"
    """
}