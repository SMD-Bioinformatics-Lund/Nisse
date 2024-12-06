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
include { GENMOD_MODELS } from './modules/genmod/genmod_models.nf'
include { GENMOD_SCORE } from './modules/genmod/genmod_score.nf'
include { GENMOD_COMPOUND } from './modules/genmod/genmod_compound.nf'
include { GENMOD_SORT } from './modules/genmod/genmod_sort.nf'

// Postprocessing
include { FILTER_VARIANTS_ON_SCORE } from './modules/postprocessing/filter_variants_on_score.nf'
include { PARSE_TOMTE_QC } from './modules/postprocessing/parse_tomte_qc.nf'
include { MAKE_SCOUT_YAML } from './modules/postprocessing/make_scout_yaml.nf'

include { validateAllParams } from './modules/utils'
include { softwareVersionsToYAML } from './modules/utils'

workflow {

    // To ponder: Do we want to validate input parameters?
    // Nice with early exit, not nice having to maintain in parallel with config
    // validateAllParams()
    ch_versions = Channel.empty()

    Channel
        .fromPath(params.csv)
        .splitCsv(header: true)
        .set { meta_ch }

    cadd_ch = Channel.of(
        tuple(
            file(params.cadd),
            file(params.cadd_tbi)
        )
    )

    vcf_ch = meta_ch.map { meta -> 
        def sample_id = meta.sample
        def variant_calls = "${params.tomte_results}/call_variants/${sample_id}_split_rmdup_info.vcf.gz"
        def variant_calls_tbi = "${variant_calls}.tbi"
        tuple(meta, file(variant_calls), file(variant_calls_tbi)) 
    }

    multiqc_ch = meta_ch.map { meta -> 
        // FIXME: Use real Tomte results when a up-to-date HG002-MultiQC run is read
        // def multiqc_summary = "${params.tomte_results}/multiqc/multiqc_data/multiqc_general_stats.txt"
        // def picard_coverage = "${params.tomte_results}/multiqc/multiqc_data/picard_rna_coverage.txt"
        def multiqc_summary = "${params.multiqc_temp}"
        def picard_coverage = "${params.picard_rna_coverage_temp}"
        tuple(meta, file(multiqc_summary), file(picard_coverage))
    }

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

    PREPROCESS(fraser_results_ch, outrider_results_ch, vcf_ch, params.hgnc_map, params.stat_col, params.stat_cutoff)
    CREATE_PED(meta_ch)

    SNV_ANNOTATE(PREPROCESS.out.vcf, cadd_ch, params.vep)
    ch_versions.mix(SNV_ANNOTATE.out.versions)
    
    SNV_SCORE(SNV_ANNOTATE.out.vcf, CREATE_PED.out.ped, params.score_config)
    ch_versions.mix(SNV_SCORE.out.versions)

    drop_results = PREPROCESS.out.fraser.join(PREPROCESS.out.outrider)
    POSTPROCESS(SNV_SCORE.out.vcf, drop_results, multiqc_ch)
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'tomte_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
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
    ch_cadd_db
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
    ADD_CADD_SCORES_TO_VCF(ch_cadd_vcf, ch_cadd_db)

    ch_versions = Channel.empty()
    ch_versions.mix(ANNOTATE_VEP.out.versions)
    ch_versions.mix(VCF_ANNO.out.versions)
    ch_versions.mix(EXTRACT_INDELS_FOR_CADD.out.versions)
    ch_versions.mix(INDEL_VEP.out.versions)
    ch_versions.mix(CALCULATE_INDEL_CADD.out.versions)
    ch_versions.mix(BGZIP_INDEL_CADD.out.versions)
    ch_versions.mix(ADD_CADD_SCORES_TO_VCF.out.versions)

    emit:
    vcf = ADD_CADD_SCORES_TO_VCF.out.vcf
    versions = ch_versions
}

workflow SNV_SCORE {
    take:
    ch_annotated_vcf
    ch_ped
    val_score_config

    main:
    GENMOD_MODELS(ch_annotated_vcf, ch_ped)
    GENMOD_SCORE(GENMOD_MODELS.out.vcf, ch_ped, val_score_config)
    GENMOD_COMPOUND(GENMOD_SCORE.out.vcf)
    GENMOD_SORT(GENMOD_COMPOUND.out.vcf)
    VCF_COMPLETION(GENMOD_SORT.out.vcf)

    ch_versions = Channel.empty()
    ch_versions.mix(GENMOD_MODELS.out.versions)
    ch_versions.mix(GENMOD_COMPOUND.out.versions)
    ch_versions.mix(GENMOD_SCORE.out.versions)
    ch_versions.mix(GENMOD_SORT.out.versions)
    ch_versions.mix(VCF_COMPLETION.out.versions)

    emit:
    vcf = VCF_COMPLETION.out.vcf
    versions = ch_versions
}

workflow POSTPROCESS {
    take:
    ch_scored_vcf
    ch_drop_results
    ch_multiqc

    main:
    FILTER_VARIANTS_ON_SCORE(ch_scored_vcf, params.score_threshold)
    MAKE_SCOUT_YAML(ch_drop_results, params.tomte_results, params.template_yaml)
    PARSE_TOMTE_QC(ch_multiqc)
}

