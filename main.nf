// A bunch of these ones will be present in fixed locations in the Tomte output

// snv_calls:
// 

def containers = ['genmod', 'vep', 'ol_wgs']
def vep_params = [
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

containers.each { param ->
    if (!params.containers.containsKey(param)) {
        params.containers[param] = null 
    }
}
vep_params.each { param ->
    if (!params.containsKey(param)) {
        params[param] = null 
    }
}
otherParams.each { param ->
    if (!params.containsKey(param)) {
        params[param] = null 
    }
}

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


def validateParams(targetParams, search_scope, type) {
    def missingParams = targetParams.findAll { !search_scope[it] }
    if (!missingParams.isEmpty()) {
        def missingList = missingParams.collect { "--${it}" }.join(", ")
        error "Error: Missing required parameter(s) in $type: ${missingList}"
    }
}

validateParams(otherParams, params, "base")
validateParams(containers, params.containers, "containers")
validateParams(vep_params, params.vep, "vep")

// OK, I can probably start stubbing this locally

workflow  {
    Channel
        .fromPath(params.csv)
        .splitCsv(header:true)
        .set { meta }

    // 

    preprocess(meta, params.hgnc_map)
        .set { out_ch }

    // snv_annotate(out_ch)
    //     .set { annotated_ch }

    // snv_score(out_ch)
    //     .set { annotated_ch }

    // postprocess(annotated_ch)
    //     .set { final_ch }

    // final_ch.view()
}

workflow preprocess {
    take:
        meta
        hgnc_map
    main:

        PREPARE_DROPS(meta, hgnc_map)
            .set { goodbye_ch }

    emit:
        goodbye_ch
}

// Something from /fs1/jakob/proj/240613_run_tomte/load_cases/4_rank_scores_run/prepare_annot_run.sh needed?
// Haven't I already implemented this?

workflow snv_annotate {
    take:
        unsure_ch
    main:
        hello(unsure_ch)
            .set { after_hello_ch }
        goodbye(after_hello_ch)
            .set { goodbye_ch }

        // What are the steps
        // [create_ped] create_ped.pl
        // [extract_indels_for_cadd] bcftools view
        // [indel_vep] vep + filter_indels.pl
        // [calculate_indel_cadd] /CADD-scripts/CADD.sh
        // [annotate_vep] vep
        // [vcfanno] vcfanno_linux64
        // [modify_vcf] modify_vcf_scout.pl
        // [mark_splice] /opt/bin/mark_spliceindels.pl
        // [add_cadd_scores_to_vcf] genmod annotate --cadd-file
        // [vcf_completion] sed + bgzip + tabix
    emit:
        goodbye_ch
}

workflow snv_score {

    // [inher_models] genmod models
    // [genmodscore] genmod score + genmod compound + sed + genmod sort

    take:
        unsure_ch
    
    main:
        hello(unsure_ch)
            .set { after_hello_ch }
    
    emit:
        goodbye_ch
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

