process MAKE_SCOUT_YAML {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.base}"

    input:
    tuple val(meta), path(fraser), path(outrider), path(vcf)
    val(tomte_results_dir)
    path(template_yaml)
    val(output_dir)

    output:
    path "${meta.sample}.yaml", emit: yaml

    script:
    """
    bash produce_yaml.sh \
        "${template_yaml}" \
        "${meta.sample}" \
        "${tomte_results_dir}" \
        "${output_dir}/${fraser}" \
        "${output_dir}/${outrider}" \
        "${vcf}" \
        "${meta.sample}.yaml"
    """

    stub:
    """
    touch "${meta.sample}_out.txt"
    """
}