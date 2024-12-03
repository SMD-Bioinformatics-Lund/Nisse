process CREATE_PED {
    tag "$meta.id"
    time '20m'
    memory '1 GB'

    input:
        tuple val(meta)

    output:
        tuple val(meta), file("${meta.id}_base.ped")

    script:
        father = "0"
        mother = "0"
        type_fa = "fa"
        type_ma = "ma"
        """
        create_ped.pl --mother $mother --father $father --group ${meta.group} --id ${meta.id} --sex ${meta.sex}
        """

    stub:
        type_fa = "fa"
        type_ma = "ma"
        """
        touch "${meta.group}_base.ped"
        touch "${meta.group}_ma.ped"
        touch "${meta.group}_fa.ped"

        echo $type_fa $type_ma > type.val
        """
}
