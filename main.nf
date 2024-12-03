include { hello } from './modules/hello'
include { goodbye } from './modules/goodbye'

include { PREPARE_DROP } from './modules/prepare_drop'
include { scout_yaml } from './modules/scout_yaml'

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
        error "Error: Missing required parameter(s) in $type: ${missingList}"
    }
}


def validateAllParams() {
    def containers = ['genmod', 'vep', 'ol_wgs']
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

workflow  {

    validateAllParams()

    Channel
        .fromPath(params.csv)
        .splitCsv(header:true)
        .set { meta_ch }

    variant_ch = meta_ch.map { meta -> tuple(meta, params.variant_calls, params.variant_calls_tbi) }

    preprocess_drop(meta_ch, params.fraser_results, params.hgnc_map)
        .set { fraser_out_ch }

    preprocess_drop(meta_ch, params.outrider_results, params.hgnc_map)
        .set { fraser_out_ch }


    snv_annotate(variant_ch)

    snv_score(variant_ch.out.ped, variant_ch.out.vcf)
}

workflow preprocess_drop {
    take:
        meta
        drop_results
        hgnc_map
    main:
        PREPARE_DROP(meta, drop_results, hgnc_map)
            .set { goodbye_ch }

    emit:
        goodbye_ch
}

workflow snv_annotate {
    take:
        ch_vcf  // channel: [mandatory] [ val(meta), path(vcf), path(vcf_tbi) ]
        ch_cadd // channel: [mandatory] [ path(cadd), path(cadd_tbi) ]

    main:

        CREATE_PED(ch_vcf[0])
            .set { ch_ped }

        // CADD indels
        EXTRACT_INDELS_FOR_CADD(ch_vcf)
        INDEL_VEP(ch_vcf).set { ch_vep_indels_only }
        CALCULATE_INDEL_CADD(ch_vep_indels_only).set { ch_cadd_indels }

        ANNOTATE_VEP(ch_vcf).set { ch_vep }
        VCF_ANNO(ch_vep).set { ch_vcf_anno }
        MODIFY_VCF(ch_vcf_anno).set { ch_scout_modified }
        MARK_SPLICE(ch_scout_modified).set { ch_mark_splice }

        ADD_CADD_SCORES_TO_VCF(ch_mark_splice, ch_cadd).set { ch_vcf_with_cadd }
        VCF_COMPLETION(ch_vcf_with_cadd).set { ch_vcf_completed }
    
    emit:
        ped = ch_ped
        vcf = ch_vcf
}

workflow snv_score {

    take:
        ch_annotated_vcf // channel: [mandatory] [ val(meta), path(vcf), path(vcf_tbi) ]
        ch_ped // channel: [mandatory] [ path(ped) ]
    
    main:
        GENMOD_MODELS(ch_annotated_vcf, ch_ped).set { ch_genmod_models }
        GENMOD_ANNOTATE(ch_genmod_models).set { ch_genmod_annotate }
        GENMOD_COMPOUND(ch_genmod_annotate).set { ch_genmod_compound }
        GENMOD_SCORE(ch_genmod_compound).set { ch_genmod_score }
    
    emit:
        ch_genmod_score
}

workflow postprocess {
    take:
        scored_vcf_ch
        csv_ch
        multiqc_ch // Both general stats and the picard
    
    main:
        filter_variants_on_score(scored_vcf_ch, params.score_threshold)

        scout_yaml(csv_ch)
            .set { after_hello_ch }

        parse_tomte_qc(multiqc_ch)
        
    emit:
        after_hello_ch
}

