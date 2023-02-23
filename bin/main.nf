#!/usr/bin/env nextflow
nextflow.preview.dsl=2

def helpMessage() {
    log.info"""
overwrite
    """.stripIndent()
}

// show help message
params.help = false
// The params scope allows you to define parameters that will be accessible in the pipeline script
if (params.help){
    helpMessage()
    exit 0
}

// Import modules from modules files
include { Overwrite_wrapper } from './modules.nf'



// Define input channels 
Star_index_Ch = Channel
            .fromPath("${params.INDEX}/star_host/")


workflow{
    Overwrite_wrapper( 
        file("${params.INPUT}"),
        file("${projectDir}/File_Compatible.py"),
        file("${projectDir}/OverwriTE_V5.sh"),
        file("${projectDir}/OverwriTE_V11.py")


    )
}