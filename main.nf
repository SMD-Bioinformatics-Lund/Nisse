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
include { BGZIP_TABIX as BGZIP_TABIX_VCF } from './modules/postprocessing/bgzip_tabix.nf'
include { BGZIP_TABIX as BGZIP_TABIX_JUNCTION_BED } from './modules/postprocessing/bgzip_tabix.nf'
include { OUTPUT_VERSIONS } from './modules/postprocessing/output_versions.nf'

include { TOMTE } from './tomte/workflows/tomte.nf'
include { IDSNP_CALL } from './modules/defined_calls.nf'
include { IDSNP_VCF_TO_JSON } from './modules/defined_calls.nf'
include { PERC_HETEROZYGOTES } from './modules/defined_calls.nf'
include { versions } from './modules/postprocessing/bgzip_tabix.nf'

workflow {

    startupMessage()

    ch_versions = Channel.empty()
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { meta ->
            // Needed for Tomte's internal workings
            // FIXME: We need to calculate the number of fq-pairs here when running many samples
            meta = meta + [ fq_pairs: 1, single_end: false, is_fastq: true, id: meta.sample ]
            meta
        }
        .set { ch_meta }

    ch_fasta_fai = Channel.of(tuple(file(params.genome_fasta), file(params.genome_fasta_fai)))

    // Either execute Tomte as part of Nisse, or start with its results folder
    if (params.run_tomte) {
        ch_meta
            .map { meta ->
                // Also needed for Tomte's internal workings
                def fastq_fw = meta.fastq_1
                def fastq_rv = meta.fastq_2
                tuple(meta, [fastq_fw, fastq_rv])
            }
            .set { ch_tomte }

        TOMTE(ch_tomte)
        ch_versions.mix(TOMTE.out.versions)

        // Tomte adds "id: meta.sample" in one step making the Nisse and Tomte
        // meta objects different
        ch_tomte_meta = TOMTE.out.bam_bai.map { it -> it[0] }

        ch_multiqc = ch_tomte_meta.combine(TOMTE.out.multiqc_data).map { meta, multiqc_folder ->
            def multiqc_summary = file("${multiqc_folder}/multiqc_general_stats.txt")
            def star_qc = file("${multiqc_folder}/multiqc_star.txt")
            def picard_coverage = file("${multiqc_folder}/picard_rna_coverage.txt")
            tuple(meta, multiqc_summary, star_qc, picard_coverage)
        }

        ch_hb_estimates = TOMTE.out.hb_estimates

        // FIXME: Is this one published to the Tomte results? No need to publish it in Nisse as well or?
        ch_junction_bed_tbi = TOMTE.out.junction_bed


        ch_ped = TOMTE.out.ped
        ch_vcf_tbi = TOMTE.out.vcf_tbi
        ch_drop_ae_out_research = TOMTE.out.drop_ae_out_research
        ch_drop_as_out_research = TOMTE.out.drop_as_out_research
        ch_bam_bai = TOMTE.out.bam_bai

    } else {
        // Creating a channel for Hb percentage from Tomte results
        ch_hb_estimates = ch_meta.map { meta ->
            def sample_id = meta.sample
            def hb_estimate_json = String.format(params.tomte_results_paths.hb_estimate, params.tomte_results, sample_id)
            tuple(meta, file(hb_estimate_json))
        }

        ch_multiqc = ch_meta.map { meta ->
            def multiqc_summary = String.format(params.tomte_results_paths.multiqc_summary, params.tomte_results)
            def star_qc = String.format(params.tomte_results_paths.star_qc, params.tomte_results)
            def picard_coverage = String.format(params.tomte_results_paths.picard_coverage, params.tomte_results)
            tuple(meta, file(multiqc_summary), file(star_qc), file(picard_coverage))
        }

        ch_vcf = ch_meta.map { meta ->
            def sample_id = meta.sample
            def variant_calls = String.format(params.tomte_results_paths.variant_calls, params.tomte_results, sample_id)
            def variant_calls_tbi = "${variant_calls}.tbi"
            tuple(meta, file(variant_calls), file(variant_calls_tbi))
        }

        ch_bam_bai = ch_meta.map { meta ->
            def sample_id = meta.sample
            def bam = String.format(params.tomte_results_paths.bam, params.tomte_results, sample_id)
            def bai = String.format(params.tomte_results_paths.bam, params.tomte_results, sample_id)
            tuple(meta, file(bam), file(bai))
        }

        if (!params.qc_only) {
            BGZIP_TABIX_VCF(ch_vcf)
            versions.mix(BGZIP_TABIX_VCF.out.versions)
            ch_vcf_tbi = BGZIP_TABIX_VCF.out.vcf_tbi

            ch_junction_bed = ch_meta.map { meta ->
                def sample_id = meta.sample
                def junction_bed = String.format(params.tomte_results_paths.junction_bed, params.tomte_results, sample_id)
                tuple(meta, file(junction_bed))
            }
            BGZIP_TABIX_JUNCTION_BED(ch_junction_bed)
            versions.mix(BGZIP_TABIX_JUNCTION_BED.out.versions)
            ch_junction_bed_tbi = BGZIP_TABIX_JUNCTION_BED.out.bed_tbi

            ch_drop_ae_out_research = ch_meta.map { meta ->
                def case_id = meta.case
                def outrider_results = String.format(params.tomte_results_paths.outrider_tsv, params.tomte_results, case_id)
                tuple(meta, file(outrider_results))
            }

            ch_drop_as_out_research = ch_meta.map { meta ->
                def case_id = meta.case
                def fraser_results = String.format(params.tomte_results_paths.fraser_tsv, params.tomte_results, case_id)
                tuple(meta, file(fraser_results))
            }
        }
    }

    NISSE_QC(
        ch_versions,
        ch_multiqc,
        ch_hb_estimates,
        ch_bam_bai,
        ch_fasta_fai,
        params.idsnps,
        params.het_calls
    )
    ch_versions = ch_versions.mix(NISSE_QC.out.versions)

    if (!params.qc_only) {
        NISSE(
            ch_versions,
            ch_meta,
            ch_junction_bed_tbi,
            ch_ped,
            ch_vcf_tbi,
            ch_drop_ae_out_research,
            ch_drop_as_out_research
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
    ch_fasta_fai
    val_idsnp_params
    val_hetcalls_params

    main:
    PERC_HETEROZYGOTES(ch_bam_bai, ch_fasta_fai, val_hetcalls_params)
    IDSNP_CALL(ch_bam_bai, ch_fasta_fai, val_idsnp_params)
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
    ch_meta
    ch_tomte_junction_bed_tbi
    ch_combined_ped
    ch_tomte_vcf_tbi
    ch_tomte_drop_ae_out_research
    ch_tomte_drop_as_out_research

    main:
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

    ch_drop_ae_per_sample = ch_meta.combine(ch_tomte_drop_as_out_research)
    ch_drop_as_per_sample = ch_meta.combine(ch_tomte_drop_ae_out_research)

    PREPROCESS(ch_drop_ae_per_sample, ch_drop_as_per_sample, ch_tomte_vcf_tbi, params.hgnc_map, params.stat_col, params.stat_cutoff)

    SNV_ANNOTATE(PREPROCESS.out.vcf, params.vep)
    ch_versions = ch_versions.mix(SNV_ANNOTATE.out.versions)

    SNV_SCORE(ch_meta, SNV_ANNOTATE.out.vcf, ch_combined_ped, params.score_config, params.score_threshold)
    ch_versions = ch_versions.mix(SNV_SCORE.out.versions)

    ch_drop_results = PREPROCESS.out.fraser.join(PREPROCESS.out.outrider)

    ch_all_result_files = ch_drop_results
        .join(SNV_SCORE.out.vcf_tbi)
        .join(ch_tomte_junction_bed_tbi)
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
    ch_meta
    ch_annotated_vcf
    ch_combined_ped
    val_score_config
    val_score_threshold

    main:
    ch_make_case_ped = ch_meta.join(ch_combined_ped)
    MAKE_CASE_PED(ch_make_case_ped)
    GENMOD_MODELS(ch_annotated_vcf, MAKE_CASE_PED.out.ped)
    GENMOD_SCORE(GENMOD_MODELS.out.vcf, MAKE_CASE_PED.out.ped, val_score_config)
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
