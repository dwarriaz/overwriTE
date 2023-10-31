#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include{
    Batch;
    Query;
    Gather;
} from './modules.nf'

workflow
{
    Batch(file(params.INPUT))
    workChannel = Batch.out
    workChannel.view
    //batchChannel = workChannel | Query
    //Gather(batchChannel.collect(),'fuckItWeBall')
}