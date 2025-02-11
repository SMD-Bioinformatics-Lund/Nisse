// Preprocessing
include { PREPARE_DROP as PREPARE_DROP_FRASER } from './modules/prepare_drop'
include { PREPARE_DROP as PREPARE_DROP_OUTRIDER } from './modules/prepare_drop'

include { BGZIP_INDEL_CADD } from './modules/bgzip_indel_cadd.nf'
include { PREPARE_VCF } from './modules/prepare_vcf.nf'

// Annotations
include { ADD_CADD_SCORES_TO_VCF } from './modules/annotate/add_cadd_scores_to_vcf.nf'
include { ANNOTATE_VEP } from './modules/annotate/annotate_vep.nf'
include { CALCULATE_INDEL_CADD } from './modules/annotate/calculate_indel_cadd.nf'
include { EXTRACT_INDELS_FOR_CADD } from './modules/annotate/extract_indels_for_cadd.nf'
include { INDEL_VEP } from './modules/annotate/indel_vep.nf'
include { MARK_SPLICE } from './modules/annotate/mark_splice.nf'
include { MODIFY_VCF } from './modules/annotate/modify_vcf.nf'
include { VCF_ANNO } from './modules/annotate/vcf_anno.nf'
include { VCF_COMPLETION } from './modules/annotate/vcf_completion.nf'

// Genmod
include { GENMOD_MODELS } from './modules/genmod/genmod_models.nf'
include { GENMOD_SCORE } from './modules/genmod/genmod_score.nf'
include { GENMOD_COMPOUND } from './modules/genmod/genmod_compound.nf'
include { GENMOD_SORT } from './modules/genmod/genmod_sort.nf'

// Postprocessing
include { FILTER_VARIANTS_ON_SCORE } from './modules/postprocessing/filter_variants_on_score.nf'
include { PARSE_TOMTE_QC } from './modules/postprocessing/parse_tomte_qc.nf'
include { MAKE_SCOUT_YAML } from './modules/postprocessing/make_scout_yaml.nf'
include { BGZIP_TABIX as BGZIP_TABIX_VCF } from './modules/postprocessing/bgzip_tabix.nf'
include { BGZIP_TABIX as BGZIP_TABIX_BED } from './modules/postprocessing/bgzip_tabix.nf'

include { BGZIP_TABIX } from './modules/postprocessing/bgzip_tabix.nf'
include { OUTPUT_VERSIONS } from './modules/postprocessing/output_versions.nf'

include { TOMTE } from './tomte/workflows/tomte.nf'

workflow {


    startupMessage()

    ch_versions = Channel.empty()
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .set { ch_meta }

    // Creating a channel for Hb percentage form Tomte results
    ch_hb_estimates = ch_meta.map { meta ->
        def sample_id = meta.sample
        def hb_estimate_json = String.format(params.tomte_results_paths.hb_estimate, params.tomte_results, sample_id)
        tuple(meta, file(hb_estimate_json))
    }

    ch_multiqc = ch_meta.map { meta ->
        def multiqc_summary = String.format(params.tomte_results_paths.multiqc_summary, params.tomte_results)
        def picard_coverage = String.format(params.tomte_results_paths.picard_coverage, params.tomte_results)
        tuple(meta, file(multiqc_summary), file(picard_coverage))
    }


    NISSE_QC(ch_versions, ch_multiqc.join(ch_hb_estimates))

    ch_versions = ch_versions.mix(NISSE_QC.out.versions)
    if (!params.qc_only) {
        NISSE(ch_versions, ch_meta)
    }
    ch_versions = ch_versions.mix(NISSE_QC.out.versions)

    ch_joined_versions = ch_versions.collect { it[1] }
    OUTPUT_VERSIONS(ch_joined_versions)

    workflow.onComplete {
        log.info("Completed without errors")
    }

    workflow.onError {
        log.error("Aborted with errors")
    }
}

workflow NISSE_QC {
    take:
    ch_versions
    ch_multiqc

    main:
    PARSE_TOMTE_QC(ch_multiqc)

    emit:
    versions = ch_versions
}

workflow NISSE {
    take:
    ch_versions
    ch_meta

    main:

    ch_tomte = ch_meta.map { meta ->
        def fastq_fw = meta.fastq_1
        def fastq_rv = meta.fastq_2

        meta = meta + [ fq_pairs: 1, single_end: false, is_fastq: true ]

        tuple(meta, [fastq_fw, fastq_rv])
    }

    TOMTE(ch_tomte)

    // ch_vcf = ch_meta.map { meta ->
    //     def sample_id = meta.sample
    //     def variant_calls = String.format(params.tomte_results_paths.variant_calls, params.tomte_results, sample_id)
    //     def variant_calls_tbi = "${variant_calls}.tbi"
    //     tuple(meta, file(variant_calls), file(variant_calls_tbi))
    // }

    // ch_junction_bed = ch_meta.map { meta ->
    //     def sample_id = meta.sample
    //     def junction_bed = String.format(params.tomte_results_paths.junction_bed, params.tomte_results, sample_id)
    //     tuple(meta, file(junction_bed))
    // }

    // ch_fraser_results = ch_meta.map { meta ->
    //     def case_id = meta.case
    //     def fraser_results = String.format(params.tomte_results_paths.fraser_tsv, params.tomte_results, case_id)
    //     tuple(meta, file(fraser_results))
    // }

    // ch_outrider_results = ch_meta.map { meta ->
    //     def case_id = meta.case
    //     def outrider_results = String.format(params.tomte_results_paths.outrider_tsv, params.tomte_results, case_id)
    //     tuple(meta, file(outrider_results))
    // }

    // NOTE: These are not accessed directly - the paths are used in the Scout yaml
    ch_tomte_raw_results = ch_meta.map { meta ->
        def sample_id = meta.sample
        def cram = String.format(params.tomte_results_paths.cram, params.tomte_results, sample_id)
        def cram_crai = String.format(params.tomte_results_paths.cram_crai, params.tomte_results, sample_id)
        def bigwig = String.format(params.tomte_results_paths.bigwig, params.tomte_results, sample_id)
        def peddy_ped = String.format(params.tomte_results_paths.peddy_ped, params.tomte_results, sample_id)
        def peddy_check = String.format(params.tomte_results_paths.peddy_check, params.tomte_results, sample_id)
        def peddy_sex = String.format(params.tomte_results_paths.peddy_sex, params.tomte_results, sample_id)
        tuple(meta, file(cram), file(cram_crai), file(bigwig), file(peddy_ped), file(peddy_check), file(peddy_sex))
    }

    PREPROCESS(TOMTE.out.drop_as_out_clinical, TOMTE.out.drop_ae_out_clinical, TOMTE.out.vcf_tbi, params.hgnc_map, params.stat_col, params.stat_cutoff)
    // PREPROCESS(ch_fraser_results, ch_outrider_results, ch_vcf, params.hgnc_map, params.stat_col, params.stat_cutoff)

    SNV_ANNOTATE(PREPROCESS.out.vcf, params.vep)
    ch_versions = ch_versions.mix(SNV_ANNOTATE.out.versions)

    SNV_SCORE(SNV_ANNOTATE.out.vcf, TOMTE.out.ped, params.score_config, params.score_threshold)
    ch_versions = ch_versions.mix(SNV_SCORE.out.versions)

    ch_drop_results = PREPROCESS.out.fraser.join(PREPROCESS.out.outrider)

    BGZIP_TABIX_BED(TOMTE.out.junction_bed)
    // BGZIP_TABIX_BED(ch_junction_bed)

    ch_all_result_files = ch_drop_results
        .join(SNV_SCORE.out.vcf_tbi)
        .join(BGZIP_TABIX_BED.out.bed_tbi)
        .join(ch_tomte_raw_results)
    MAKE_SCOUT_YAML(ch_all_result_files, params.tomte_results, params.outdir, params.phenotype, params.tissue)

    emit:
    versions = ch_versions
}

workflow PREPROCESS {
    take:
    ch_fraser_results
    ch_outrider_results
    ch_vcf
    val_hgnc_map
    val_stat_col
    val_stat_cutoff

    main:
    PREPARE_DROP_FRASER("FRASER", ch_fraser_results, val_hgnc_map, val_stat_col, val_stat_cutoff)
    PREPARE_DROP_OUTRIDER("OUTRIDER", ch_outrider_results, val_hgnc_map, val_stat_col, val_stat_cutoff)
    PREPARE_VCF(ch_vcf)

    emit:
    vcf = PREPARE_VCF.out.vcf
    fraser = PREPARE_DROP_FRASER.out.drop
    outrider = PREPARE_DROP_OUTRIDER.out.drop
}

workflow SNV_ANNOTATE {
    take:
    ch_vcf
    val_vep_params

    main:
    ANNOTATE_VEP(ch_vcf, val_vep_params)
    VCF_ANNO(ANNOTATE_VEP.out.vcf, val_vep_params)
    MODIFY_VCF(VCF_ANNO.out.vcf)
    MARK_SPLICE(MODIFY_VCF.out.vcf)

    EXTRACT_INDELS_FOR_CADD(ch_vcf)
    INDEL_VEP(EXTRACT_INDELS_FOR_CADD.out.vcf, val_vep_params)
    CALCULATE_INDEL_CADD(INDEL_VEP.out.vcf)
    BGZIP_INDEL_CADD(CALCULATE_INDEL_CADD.out.vcf)

    ch_cadd_vcf = MARK_SPLICE.out.vcf.join(BGZIP_INDEL_CADD.out.cadd)
    ADD_CADD_SCORES_TO_VCF(ch_cadd_vcf)

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(ANNOTATE_VEP.out.versions)
    ch_versions = ch_versions.mix(VCF_ANNO.out.versions)
    ch_versions = ch_versions.mix(EXTRACT_INDELS_FOR_CADD.out.versions)
    ch_versions = ch_versions.mix(INDEL_VEP.out.versions)
    ch_versions = ch_versions.mix(CALCULATE_INDEL_CADD.out.versions)
    ch_versions = ch_versions.mix(BGZIP_INDEL_CADD.out.versions)
    ch_versions = ch_versions.mix(ADD_CADD_SCORES_TO_VCF.out.versions)

    emit:
    vcf = ADD_CADD_SCORES_TO_VCF.out.vcf
    versions = ch_versions
}

workflow SNV_SCORE {
    take:
    ch_annotated_vcf
    ch_ped
    val_score_config
    val_score_threshold

    main:
    GENMOD_MODELS(ch_annotated_vcf, ch_ped)
    GENMOD_SCORE(GENMOD_MODELS.out.vcf, ch_ped, val_score_config)
    GENMOD_COMPOUND(GENMOD_SCORE.out.vcf)
    GENMOD_SORT(GENMOD_COMPOUND.out.vcf)
    VCF_COMPLETION(GENMOD_SORT.out.vcf)
    FILTER_VARIANTS_ON_SCORE(VCF_COMPLETION.out.vcf, val_score_threshold)
    BGZIP_TABIX_VCF(FILTER_VARIANTS_ON_SCORE.out.vcf)

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(GENMOD_MODELS.out.versions)
    ch_versions = ch_versions.mix(GENMOD_COMPOUND.out.versions)
    ch_versions = ch_versions.mix(GENMOD_SCORE.out.versions)
    ch_versions = ch_versions.mix(GENMOD_SORT.out.versions)
    ch_versions = ch_versions.mix(VCF_COMPLETION.out.versions)
    ch_versions = ch_versions.mix(BGZIP_TABIX_VCF.out.versions)

    emit:
    vcf_tbi = BGZIP_TABIX_VCF.out.vcf_tbi
    versions = ch_versions
}

def startupMessage() {
    print("Starting Nisse")
    print("Output dir: ${params.outdir}")
}
