process PREPARE_VCF {

    tag "${meta.sample}"
    label "process_low"
    container "${params.containers.bcftools}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple path("${meta.sample}_prepared.vcf.gz"), path("${meta.sample}_prepared.vcf.gz.tbi"), emit: vcf

    script:
    """
    # We need a dummy Pathogenicity attribute to not crash without MT                                                                   ────────────────────────────────────────────────────────────────────────────
    # Version=3 in the VDB crashed genmod

    zcat ${vcf} |\
        sed "/#CHROM/i ##INFO=<ID=Pathogenicity,Number=1,Type=Integer,Description=\"Placeholder to not crash genmod\">" |\
        sed "/ID=VDB/ s/,Version=\"3\"//" |\
        grep -P "^#|^chr[1-9]\b|^chr[12][0-9]\b|^chr[MXY]\b" > ${meta.sample}_prepared.vcf
    bgzip "${meta.sample}_prepared.vcf"
    tabix "${meta.sample}_prepared.vcf.gz"
    """

    stub:
    """
    touch "${meta.sample}_prepared.vcf.gz"
    touch "${meta.sample}_prepared.vcf.gz.tbi"
    """
}