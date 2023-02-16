#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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
include { Trimming_FastP } from './illumina_modules.nf'



// Define input channels 
Star_index_Ch = Channel
            .fromPath("${params.INDEX}/star_host/")


workflow{}