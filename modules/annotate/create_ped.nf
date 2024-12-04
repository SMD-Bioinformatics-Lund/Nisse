process CREATE_PED {
    tag "${meta.sample}"
    time '20m'
    memory '1 GB'

    input:
        val(meta)

    output:
        tuple val(meta), file("${meta.case}_base.ped"), emit: ped

    script:
        father = "0"
        mother = "0"
        type_fa = "fa"
        type_ma = "ma"
        def group = "${meta.case}"
        def sample = "${meta.sample}"
        def sex = "${meta.sex}"
        """
        create_ped.pl --mother "${mother}" --father "${father}" --group "${group}" --id "${sample}" --sex "${sex}"
        """

    stub:
        type_fa = "fa"
        type_ma = "ma"
        def group = "${meta.case}"
        """
        touch "${group}_base.ped"
        touch "${group}_ma.ped"
        touch "${group}_fa.ped"

        echo $type_fa $type_ma > ${group}_base.ped
        """
}
