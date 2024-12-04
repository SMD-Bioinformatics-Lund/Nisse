process CREATE_PED {

    tag "${meta.sample}"
	label "process_low"
	container "${params.containers.base}"

    input:
        val(meta)

    output:
        tuple val(meta), file("${meta.sample}_base.ped"), emit: ped

    script:
        father = "0"
        mother = "0"
        type_fa = "fa"
        type_ma = "ma"
        """
        create_ped.pl --mother "${mother}" --father "${father}" --group "${meta.sample}" --id "${meta.sample}" --sex "${meta.sex}"
        """

    stub:
        type_fa = "fa"
        type_ma = "ma"
        """
        touch "${meta.sample}_base.ped"
        touch "${meta.sample}_ma.ped"
        touch "${meta.sample}_fa.ped"

        echo $type_fa $type_ma > ${meta.sample}_base.ped
        """
}
