// Preprocessing
include { PREPARE_DROP as PREPARE_DROP_FRASER } from './modules/prepare_drop'
include { PREPARE_DROP as PREPARE_DROP_OUTRIDER } from './modules/prepare_drop'

include { BGZIP_INDEL_CADD } from './modules/bgzip_indel_cadd.nf'
include { PREPARE_VCF } from './modules/prepare_vcf.nf'

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
    // FIXME: Early verification of the input files

    Channel
        .fromPath(params.csv)
        .splitCsv(header: true)
        .set { meta_ch }

    vcf_ch = meta_ch.map { meta -> 
        def sample_id = meta.sample
        def variant_calls = "${params.tomte_results}/call_variants/${sample_id}_split_rmdup_info.vcf.gz"
        def variant_calls_tbi = "${variant_calls}.tbi"
        tuple(meta, file(variant_calls), file(variant_calls_tbi)) 
    }

    multiqc_ch = meta_ch.map { meta -> 
        // FIXME: Restore when up to date data is generated
        // def multiqc_summary = "${params.tomte_results}/multiqc/multiqc_data/multiqc_general_stats.txt"
        // def picard_coverage = "${params.tomte_results}/multiqc/multiqc_data/picard_rna_coverage.txt"
        def multiqc_summary = "${params.multiqc_temp}"
        def picard_coverage = "${params.picard_rna_coverage_temp}"
        tuple(meta, file(multiqc_summary), file(picard_coverage))
    }

    Channel
        .fromPath(params.hgnc_map)
        .set { hgnc_map_ch }

    fraser_results_ch = meta_ch.map { meta ->
        def case_id = meta.case
        def fraser_results = "${params.tomte_results}/analyse_transcripts/drop/${case_id}_fraser_top_hits_research.tsv"
        tuple(meta, file(fraser_results))
    }

    outrider_results_ch = meta_ch.map { meta ->
        def case_id = meta.case
        def outrider_results = "${params.tomte_results}/analyse_transcripts/drop/${case_id}_outrider_top_hits_research.tsv"
        tuple(meta, file(outrider_results))
    }

    // FIXME: Look into DROP processing at the end
    PREPROCESS(fraser_results_ch, outrider_results_ch, hgnc_map_ch, vcf_ch)

    Channel
        .of(tuple(params.cadd, params.cadd_tbi))
        .set { cadd_ch }

    Channel
        .fromPath(params.score_config)
        .set { score_config_ch }

    CREATE_PED(meta_ch)
    SNV_ANNOTATE(PREPROCESS.out.vcf)
    SNV_SCORE(SNV_ANNOTATE.out.vcf, CREATE_PED.out.ped, score_config_ch)

    drop_results = PREPROCESS.out.fraser.join(PREPROCESS.out.outrider)

    POSTPROCESS(SNV_SCORE.out.vcf, drop_results, multiqc_ch)

    // FIXME: Collect versions
}

    // ch_annotated_vcf // channel: [mandatory] [ val(meta), path(vcf), path(vcf_tbi) ]
    // ch_ped // channel: [mandatory] [ path(ped) ]
    // ch_score_config // channel: [mandatory] [ path(score_config) ]

workflow PREPROCESS {
    take:
    ch_fraser_results
    ch_outrider_results
    ch_hgnc_map
    ch_vcf

    main:
    PREPARE_DROP_FRASER("FRASER", ch_fraser_results, ch_hgnc_map)
    PREPARE_DROP_OUTRIDER("OUTRIDER", ch_outrider_results, ch_hgnc_map)
    PREPARE_VCF(ch_vcf)

    emit:
    vcf = PREPARE_VCF.out.vcf
    fraser = PREPARE_DROP_FRASER.out.drop
    outrider = PREPARE_DROP_OUTRIDER.out.drop
}

workflow SNV_ANNOTATE {
    take:
    ch_vcf // channel: [mandatory] [ val(meta), path(vcf), path(vcf_tbi) ]

    main:
    ANNOTATE_VEP(ch_vcf)
    VCF_ANNO(ANNOTATE_VEP.out.vcf)
    MODIFY_VCF(VCF_ANNO.out.vcf)
    MARK_SPLICE(MODIFY_VCF.out.vcf)

    // CADD indels
    EXTRACT_INDELS_FOR_CADD(ch_vcf)
    INDEL_VEP(EXTRACT_INDELS_FOR_CADD.out.vcf)
    CALCULATE_INDEL_CADD(INDEL_VEP.out.vcf)
    BGZIP_INDEL_CADD(CALCULATE_INDEL_CADD.out.vcf)

    cadd_ch = MARK_SPLICE.out.vcf.join(BGZIP_INDEL_CADD.out.cadd)
    ADD_CADD_SCORES_TO_VCF(cadd_ch)

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

workflow POSTPROCESS {
    take:
    scored_vcf_ch
    ch_drop_results
    multiqc_ch // Both general stats and the picard

    main:
    FILTER_VARIANTS_ON_SCORE(scored_vcf_ch, params.score_threshold)
    MAKE_SCOUT_YAML(ch_drop_results, params.tomte_results, params.template_yaml)
    PARSE_TOMTE_QC(multiqc_ch)
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
