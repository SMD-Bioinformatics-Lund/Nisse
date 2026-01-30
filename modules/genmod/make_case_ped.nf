process MAKE_CASE_PED {

    tag "${meta.sample}"
    label "process_low"

    input:
    tuple val(meta)
    path(family_all_cases)

    output:
    tuple val(meta), path("${meta.case}_per_case.ped"), emit: ped

    script:
    """
    grep -E "^#|^${meta.case}\\s" "${family_all_cases}" > "${meta.case}_per_case.ped"
    """

    stub:
    """
    touch "${meta.case}_per_case.ped"
    """
}
