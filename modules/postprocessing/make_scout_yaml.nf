process MAKE_SCOUT_YAML {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.base}"

    input:
    tuple val(meta), path(fraser), path(outrider)
    path(results_dir)
    path(template_yaml)

    output:
    path "${meta.sample}.yaml", emit: yaml

    script:
    """
    bash produce_yaml.sh \
        "${template_yaml}" \
        "${meta.sample}" \
        "${results_dir}" \
        "${fraser}" \
        "${outrider}" \
        "${meta.sample}.yaml"
    """

    stub:
    """
    touch "${meta.sample}_out.txt"
    """
}