process FILTER_VARIANTS_ON_SCORE {
    input:
    tuple path(scored_vcf), path(scored_vcf_tbi)
    val(threshold)

    output:
    path 'hello_out.txt'

    script:
    """
    bash hello.sh ${scored_vcf} > "hello_out.txt"
    """

    stub:
    """
    touch hello_out.txt
    """
}