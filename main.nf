nextflow.enable.dsl=2

include { hello, goodbye } from './modules'

workflow  {
    Channel
        .fromPath(params.csv)
        .set { csv_ch }

    hello(csv_ch)
        .set { after_hello_ch }
    
    goodbye(after_hello_ch)
        .set { goodbye_ch }

    goodbye_ch.view()
    
}