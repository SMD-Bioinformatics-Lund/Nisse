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
include { MAKE_CASE_PED } from './modules/genmod/make_case_ped.nf'

// Postprocessing
include { FILTER_VARIANTS_ON_SCORE } from './modules/postprocessing/filter_variants_on_score.nf'
include { PARSE_QC_FOR_CDM } from './modules/postprocessing/parse_qc_for_cdm.nf'
include { MAKE_SCOUT_YAML } from './modules/postprocessing/make_scout_yaml.nf'
include { BGZIP as BGZIP_JUNCTION_BED } from './modules/postprocessing/bgzip.nf'
include { TABIX as TABIX_JUNCTION_BED } from './modules/postprocessing/tabix.nf'
// include { BGZIP_TABIX as BGZIP_TABIX_JUNCTION_BED } from './modules/postprocessing/bgzip_tabix.nf'
include { BGZIP_TABIX as BGZIP_TABIX_VCF } from './modules/postprocessing/bgzip_tabix.nf'
include { OUTPUT_VERSIONS } from './modules/postprocessing/output_versions.nf'

include { TOMTE } from './tomte/workflows/tomte.nf'
include { IDSNP_CALL } from './modules/defined_calls.nf'
include { IDSNP_VCF_TO_JSON } from './modules/defined_calls.nf'
include { PERC_HETEROZYGOTES } from './modules/defined_calls.nf'
include { versions } from './modules/postprocessing/bgzip_tabix.nf'
include { PIPELINE_INITIALISATION } from './tomte/subworkflows/local/utils_nfcore_tomte_pipeline/main.nf'
include { CREATE_PED } from './modules/annotate/create_ped.nf'

def join_on_sample(ch1, ch2) {
    def mapped1 = ch1.map { tuple -> [tuple[0].sample, tuple] }
    def mapped2 = ch2.map { tuple -> [tuple[0].sample, tuple] }
    return mapped1
        .join(mapped2)
        .map { key_values ->
            def _key = key_values[0]
            def ch1_values = key_values[1]
            def ch2_values = key_values[2]
            [ch1_values[0]] + ch1_values.drop(1) + ch2_values.drop(1)
        }
}

workflow {

    startupMessage(params.show_params)

    ch_versions = channel.empty()
    channel.fromPath(params.input)
        .splitCsv(header: true)
        .set { ch_meta_nisse }

    // Either execute Tomte as part of Nisse, or start with its results folder
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    TOMTE(PIPELINE_INITIALISATION.out.samplesheet)

    ch_versions.mix(TOMTE.out.versions)

    // Tomte adds "id: meta.sample" in one step making the Nisse and Tomte
    // meta objects different
    ch_tomte_meta = TOMTE.out.bam_bai.map { it -> it[0] }

    ch_multiqc = ch_tomte_meta
        .combine(TOMTE.out.multiqc_data)
        .map { meta, multiqc_folder ->
            def multiqc_summary = file("${multiqc_folder}/multiqc_general_stats.txt")
            def star_qc = file("${multiqc_folder}/multiqc_star.txt")
            def picard_coverage = file("${multiqc_folder}/picard_rna_coverage.txt")
            tuple(meta, multiqc_summary, star_qc, picard_coverage)
        }

    ch_junction_bed = ch_meta_nisse.map { meta ->
        def sample_id = meta.sample
        def junction_bed_gz = String.format(params.tomte_results_paths.junction_bed_gz, params.outdir, sample_id)
        def junction_bed_gz_tbi = String.format(params.tomte_results_paths.junction_bed_gz_tbi, params.outdir, sample_id)
        tuple(meta, file(junction_bed_gz), file(junction_bed_gz_tbi))
    }

    ch_hb_estimates = TOMTE.out.hb_estimates

    // ch_vcf = TOMTE.out.vcf_tbi
    // ch_ped_tomte = TOMTE.out.ped
    ch_vcf_tbi_tomte = TOMTE.out.vcf_tbi
    ch_drop_ae_out_research_tomte = TOMTE.out.drop_ae_out_research
    ch_drop_as_out_research_tomte = TOMTE.out.drop_as_out_research
    ch_bam_bai_tomte = TOMTE.out.bam_bai

    ch_ped_nisse = CREATE_PED(ch_meta_nisse)

    NISSE_QC(
        ch_versions,
        ch_multiqc,
        ch_hb_estimates,
        ch_bam_bai_tomte,
        params.idsnps,
        params.het_calls,
    )
    ch_versions = ch_versions.mix(NISSE_QC.out.versions)

    if (!params.qc_only) {
        NISSE(
            ch_versions,
            ch_meta_nisse,
            ch_junction_bed,
            ch_ped_nisse,
            ch_vcf_tbi_tomte,
            ch_drop_ae_out_research_tomte,
            ch_drop_as_out_research_tomte,
        )
    }
    ch_versions = ch_versions.mix(NISSE.out.versions)

    // Join the paths, skipping the meta value
    ch_joined_versions = ch_versions.collect { it -> it[1] }
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
    ch_hb_estimates
    ch_bam_bai
    val_idsnp_params
    val_hetcalls_params

    main:
    PERC_HETEROZYGOTES(ch_bam_bai, val_hetcalls_params)
    IDSNP_CALL(ch_bam_bai, val_idsnp_params)
    IDSNP_VCF_TO_JSON(IDSNP_CALL.out.vcf)

    ch_qc = ch_multiqc
    ch_qc = ch_qc.join(ch_hb_estimates)
    ch_qc = ch_qc.join(PERC_HETEROZYGOTES.out.vcf)

    PARSE_QC_FOR_CDM(ch_qc)

    ch_versions = ch_versions.mix(PERC_HETEROZYGOTES.out.versions)
    ch_versions = ch_versions.mix(IDSNP_CALL.out.versions)
    ch_versions = ch_versions.mix(IDSNP_VCF_TO_JSON.out.versions)
    ch_versions = ch_versions.mix(PARSE_QC_FOR_CDM.out.versions)

    emit:
    versions = ch_versions
}

workflow NISSE {
    take:
    ch_versions
    ch_meta_nisse
    ch_junction_bed_tbi_tomte
    ch_ped_nisse
    ch_tomte_vcf_tbi_tomte
    ch_tomte_drop_ae_out_research_tomte
    ch_tomte_drop_as_out_research_tomte

    main:
    // NOTE: These are not accessed directly - the paths are used in the Scout yaml
    ch_tomte_raw_results = ch_meta_nisse.map { meta ->
        def sample_id = meta.sample
        def cram = String.format(params.tomte_results_paths.cram, params.outdir, sample_id)
        def cram_crai = String.format(params.tomte_results_paths.cram_crai, params.outdir, sample_id)
        def bigwig = String.format(params.tomte_results_paths.bigwig, params.outdir, sample_id)
        def peddy_ped = String.format(params.tomte_results_paths.peddy_ped, params.outdir, sample_id)
        def peddy_check = String.format(params.tomte_results_paths.peddy_check, params.outdir, sample_id)
        def peddy_sex = String.format(params.tomte_results_paths.peddy_sex, params.outdir, sample_id)
        tuple(meta, file(cram), file(cram_crai), file(bigwig), file(peddy_ped), file(peddy_check), file(peddy_sex))
    }

    ch_drop_ae_per_sample = ch_meta_nisse.combine(ch_tomte_drop_as_out_research_tomte)
    ch_drop_as_per_sample = ch_meta_nisse.combine(ch_tomte_drop_ae_out_research_tomte)

    PREPROCESS(ch_drop_ae_per_sample, ch_drop_as_per_sample, ch_tomte_vcf_tbi_tomte, params.hgnc_map, params.stat_col, params.stat_cutoff)

    SNV_ANNOTATE(PREPROCESS.out.vcf, params.vep)
    ch_versions = ch_versions.mix(SNV_ANNOTATE.out.versions)

    SNV_SCORE(ch_meta_nisse, SNV_ANNOTATE.out.vcf, ch_ped_nisse, params.score_config, params.score_threshold)
    ch_versions = ch_versions.mix(SNV_SCORE.out.versions)

    ch_drop_results = PREPROCESS.out.fraser.join(PREPROCESS.out.outrider)

    // Regular join did not work due to differing meta objects
    ch_2 = join_on_sample(ch_drop_results, SNV_SCORE.out.vcf_tbi)
    ch_3 = join_on_sample(ch_2, ch_junction_bed_tbi_tomte)
    ch_all_result_files = join_on_sample(ch_3, ch_tomte_raw_results)

    MAKE_SCOUT_YAML(ch_all_result_files, params.outdir, params.nisse_outdir, params.phenotype, params.tissue)

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

    ch_versions = channel.empty()
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
    ch_meta_nisse
    ch_annotated_vcf
    ch_ped_nisse
    val_score_config
    val_score_threshold

    main:
    MAKE_CASE_PED(ch_meta_nisse, ch_ped_nisse)

    ch_annotated_vcf_ped = ch_annotated_vcf.join(MAKE_CASE_PED.out.ped)

    GENMOD_MODELS(ch_annotated_vcf_ped)

    GENMOD_SCORE(GENMOD_MODELS.out.vcf_ped, val_score_config)
    GENMOD_COMPOUND(GENMOD_SCORE.out.vcf)
    GENMOD_SORT(GENMOD_COMPOUND.out.vcf)
    VCF_COMPLETION(GENMOD_SORT.out.vcf)
    FILTER_VARIANTS_ON_SCORE(VCF_COMPLETION.out.vcf, val_score_threshold)
    BGZIP_TABIX_VCF(FILTER_VARIANTS_ON_SCORE.out.vcf)

    ch_versions = channel.empty()
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

def startupMessage(showParams) {
    print("Starting Nisse")
    print("Output dir: ${params.outdir}")

    if (showParams) {
        def prettyParams = params.sort().collect { k, v -> "${k}: ${v}" }.join('\n')
        log.info("Workflow params:\n${prettyParams}")
    }
}
