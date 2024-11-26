nextflow.enable.dsl=2

def requiredParams = ['csv', 'annot', 'score_thres']
requiredParams.each { param ->
    if (!params.containsKey(param)) {
        params[param] = null 
    }
}

include { hello } from './modules/hello'
include { goodbye } from './modules/goodbye'

include { prepare_drop } from './modules/prepare_drop'
include { scout_yaml } from './modules/scout_yaml'

def validateParams(requiredParams) {
    def missingParams = requiredParams.findAll { !params[it] }
    if (!missingParams.isEmpty()) {
        def missingList = missingParams.collect { "--${it}" }.join(", ")
        error "Error: Missing required parameter(s): ${missingList}"
    }
}

validateParams(requiredParams)

// OK, I can probably start stubbing this locally

workflow  {
    Channel
        .fromPath(params.csv)
        .set { csv_ch }

    Channel
        .fromPath(params.annot)
        .set { annot_ch }

    preprocess(csv_ch)
        .set { out_ch }

    annotate(out_ch)
        .set { annotated_ch }

    postprocess(annotated_ch)
        .set { final_ch }

    final_ch.view()
    
}

workflow preprocess {
    take:
        unsure_ch
    
    main:

        prepare_drop(unsure_ch)
            .set { goodbye_ch }

    emit:
        goodbye_ch
}

workflow annotate {
    take:
        unsure_ch
    
    main:
        hello(unsure_ch)
            .set { after_hello_ch }
        
        goodbye(after_hello_ch)
            .set { goodbye_ch }        

    emit:
        goodbye_ch
}

workflow postprocess {
    take:
        scored_vcf_ch
        csv_ch
    
    // Filter out variants with higher scores
    // Produce Scout yaml
    // parse_tomte.py

    main:
        filter_variants_on_score(scored_vcf_ch, params.score_threshold)

        scout_yaml(csv_ch)
            .set { after_hello_ch }
        
    emit:
        after_hello_ch
}
