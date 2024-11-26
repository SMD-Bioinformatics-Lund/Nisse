process CREATE_PED {
    tag "$meta.id"
    time '20m'
    publishDir "${OUTDIR}/ped", mode: 'copy' , overwrite: 'true'
    memory '1 GB'

    input:
        set meta
    
    output:
        set meta, file("${meta.id}_base.ped")

    script:
        father = "0"
        mother = "0"
        type_fa = "fa"
        type_ma = "ma"
        """
        create_ped.pl --mother $mother --father $father --group $group --id $id --sex $sex
        """

    stub:
        type_fa = "fa"
        type_ma = "ma"
        """
        touch "${group}_base.ped"
        touch "${group}_ma.ped"
        touch "${group}_fa.ped"

        echo $type_fa $type_ma > type.val
        """
}
