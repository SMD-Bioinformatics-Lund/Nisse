nextflow.enable.dsl=2

params.csv = null
params.annot = null

include { hello } from './modules/hello'
include { goodbye } from './modules/goodbye'

def validateParams(requiredParams) {
    def missingParams = requiredParams.findAll { !params[it] }
    if (!missingParams.isEmpty()) {
        def missingList = missingParams.collect { "--${it}" }.join(", ")
        error "Error: Missing required parameter(s): ${missingList}"
    }
}

def requiredParams = ['csv', 'annot']
validateParams(requiredParams)

workflow  {
    Channel
        .fromPath(params.csv)
        .set { csv_ch }

    Channel
        .fromPath(params.annot)
        .set { annot_ch }

    // Verify input parameters

    // prepare_drop.sh
    // Annotation run
    // Filter out variants with higher scores
    // Produce Scout yaml
    // parse_tomte.py

    hello_goodbye(csv_ch)
        .set { goodbye_ch }

    // hello(csv_ch)
    //     .set { after_hello_ch }
    
    // goodbye(after_hello_ch)
    //     .set { goodbye_ch }

    goodbye_ch.view()
    
}

workflow hello_goodbye {
    take:
        csv_ch
    
    main:
        hello(csv_ch)
            .set { after_hello_ch }
        
        goodbye(after_hello_ch)
            .set { goodbye_ch }
        
    emit:
        goodbye_ch
}
