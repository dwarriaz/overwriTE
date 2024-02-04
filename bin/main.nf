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
    Batch.out.flatten().view() 
    
    /*Query(Batch.out.flatten())
    Gather(Query.out.collect(),'overwriTE.table.csv')
    */
}

