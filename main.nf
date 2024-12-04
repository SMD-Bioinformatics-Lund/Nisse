// Preprocessing
include { PREPARE_DROP as PREPARE_DROP_FRASER } from './modules/prepare_drop'
include { PREPARE_DROP as PREPARE_DROP_OUTRIDER } from './modules/prepare_drop'

// Annotations
include { ADD_CADD_SCORES_TO_VCF } from './modules/annotate/add_cadd_scores_to_vcf.nf'
include { ANNOTATE_VEP } from './modules/annotate/annotate_vep.nf'
include { CALCULATE_INDEL_CADD } from './modules/annotate/calculate_indel_cadd.nf'
include { CREATE_PED } from './modules/annotate/create_ped.nf'
include { EXTRACT_INDELS_FOR_CADD } from './modules/annotate/extract_indels_for_cadd.nf'
include { INDEL_VEP } from './modules/annotate/indel_vep.nf'
include { MARK_SPLICE } from './modules/annotate/mark_splice.nf'
include { MODIFY_VCF } from './modules/annotate/modify_vcf.nf'
include { VCF_ANNO } from './modules/annotate/vcf_anno.nf'
include { VCF_COMPLETION } from './modules/annotate/vcf_completion.nf'

// Genmod
include { GENMOD_ANNOTATE } from './modules/genmod/genmod_annotate.nf'
include { GENMOD_MODELS } from './modules/genmod/genmod_models.nf'
include { GENMOD_SCORE } from './modules/genmod/genmod_score.nf'
include { GENMOD_COMPOUND } from './modules/genmod/genmod_compound.nf'

// Postprocessing
include { FILTER_VARIANTS_ON_SCORE } from './modules/postprocessing/filter_variants_on_score.nf'
include { PARSE_TOMTE_QC } from './modules/postprocessing/parse_tomte_qc.nf'
include { MAKE_SCOUT_YAML } from './modules/postprocessing/make_scout_yaml.nf'

workflow {

    validateAllParams()

    // FIXME: Check that the input CSV has only one line
    fraser_results = "${params.tomte_results}/analyse_transcripts/drop/${params.case_id}_fraser_top_hits_research.tsv"
    outrider_results = "${params.tomte_results}/analyse_transcripts/drop/${params.case_id}_outrider_top_hits_research.tsv"
    variant_calls = "${params.tomte_results}/call_variants/${params.sample_id}_split_rmdup_info.vcf.gz"
    variant_calls_tbi = "${variant_calls}.tbi"


    Channel
        .fromPath(params.csv)
        .splitCsv(header: true)
        .set { meta_ch }

    vcf_ch = meta_ch.map { meta -> tuple(meta, variant_calls, variant_calls_tbi) }

    Channel
        .fromPath(params.hgnc_map)
        .set { hgnc_map_ch }

    Channel
        .fromPath(fraser_results)
        .set { fraser_results_ch }

    Channel
        .fromPath(outrider_results)
        .set { outrider_results_ch }

    // FIXME: Look into DROP processing at the end
    // preprocess(meta_ch, fraser_results_ch, outrider_results_ch, hgnc_map_ch)

    Channel
        .of(tuple(params.cadd, params.cadd_tbi))
        .set { cadd_ch }

    Channel
        .fromPath(params.score_config)
        .set { score_config_ch }

    CREATE_PED(meta_ch)
    SNV_ANNOTATE(vcf_ch, cadd_ch)
    SNV_SCORE(SNV_ANNOTATE.out.vcf, CREATE_PED.out.ped, score_config_ch)
}

    // ch_annotated_vcf // channel: [mandatory] [ val(meta), path(vcf), path(vcf_tbi) ]
    // ch_ped // channel: [mandatory] [ path(ped) ]
    // ch_score_config // channel: [mandatory] [ path(score_config) ]

workflow preprocess {
    take:
    ch_meta
    ch_fraser_results
    ch_outrider_results
    ch_hgnc_map

    main:
    PREPARE_DROP_FRASER(ch_meta, "FRASER", ch_fraser_results, ch_hgnc_map).set { fraser_ch }

    PREPARE_DROP_OUTRIDER(ch_meta, "OUTRIDER", ch_outrider_results, ch_hgnc_map).set { outrider_ch }

    emit:
    fraser_ch
    outrider_ch
}

workflow SNV_ANNOTATE {
    take:
    ch_vcf // channel: [mandatory] [ val(meta), path(vcf), path(vcf_tbi) ]
    ch_cadd // channel: [mandatory] [ path(cadd), path(cadd_tbi) ]

    main:

    // CADD indels
    EXTRACT_INDELS_FOR_CADD(ch_vcf)
    INDEL_VEP(EXTRACT_INDELS_FOR_CADD.out.vcf)
    CALCULATE_INDEL_CADD(INDEL_VEP.out.vcf)

    ANNOTATE_VEP(CALCULATE_INDEL_CADD.out.vcf)
    VCF_ANNO(ANNOTATE_VEP.out.vcf)
    MODIFY_VCF(VCF_ANNO.out.vcf)
    MARK_SPLICE(MODIFY_VCF.out.vcf)

    ADD_CADD_SCORES_TO_VCF(MARK_SPLICE.out.vcf, ch_cadd)

    emit:
    vcf = ADD_CADD_SCORES_TO_VCF.out.vcf
}

workflow SNV_SCORE {
    take:
    ch_annotated_vcf // channel: [mandatory] [ val(meta), path(vcf), path(vcf_tbi) ]
    ch_ped // channel: [mandatory] [ path(ped) ]
    ch_score_config // channel: [mandatory] [ path(score_config) ]


    main:

    GENMOD_MODELS(ch_annotated_vcf, ch_ped)
    GENMOD_ANNOTATE(GENMOD_MODELS.out.vcf)
    GENMOD_COMPOUND(GENMOD_ANNOTATE.out.vcf)

    GENMOD_SCORE(GENMOD_COMPOUND.out.vcf, ch_ped, ch_score_config)
    VCF_COMPLETION(GENMOD_SCORE.out.vcf)

    emit:
    vcf = VCF_COMPLETION.out.vcf
}

workflow postprocess {
    take:
    scored_vcf_ch
    csv_ch
    multiqc_ch // Both general stats and the picard

    main:
    FILTER_VARIANTS_ON_SCORE(scored_vcf_ch, params.score_threshold)

    MAKE_SCOUT_YAML(csv_ch).set { after_hello_ch }

    PARSE_TOMTE_QC(multiqc_ch)

    emit:
    after_hello_ch
}

// OK some thinking
// In point is output from Tomte
// 1. SNV calls on RNA-seq
// 2. DROP results
// In reality the DROP results will be for a single sample, isn't it?
// We can maybe assume that pre-processing here

// OK, and now I can start with drafting the stub run


def assignDefaultParams(target_params, user_params) {
    target_params.each { param ->
        if (!user_params.containers.containsKey(param)) {
            user_params.containers[param] = null
        }
    }
}

def validateParams(targetParams, search_scope, type) {
    def missingParams = targetParams.findAll { !search_scope[it] }
    if (!missingParams.isEmpty()) {
        def missingList = missingParams.collect { "--${it}" }.join(", ")
        error("Error: Missing required parameter(s) in ${type}: ${missingList}")
    }
}


def validateAllParams() {
    def containers = ['genmod', 'vep', 'cadd', 'base']
    def vepParams = [
        'VEP_SYNONYMS',
        'VEP_FASTA',
        'VEP_CACHE',
        'VEP_PLUGINS',
        'VEP_TRANSCRIPT_DISTANCE',
        'CADD',
        'MAXENTSCAN',
        'DBNSFP',
        'GNOMAD_EXOMES',
        'GNOMAD_GENOMES',
        'GNOMAD_MT',
        'PHYLOP',
        'PHASTCONS'
    ]

    def otherParams = ['csv', 'score_thres', 'snv_calls']

    assignDefaultParams(containers, params)
    assignDefaultParams(vepParams, params)
    assignDefaultParams(otherParams, params)

    validateParams(otherParams, params, "base")
    validateParams(containers, params.containers, "containers")
    validateParams(vepParams, params.vep, "vep")
}
