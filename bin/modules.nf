process Batch{
    input: 
        path inputFile
    
    output: 
        path "batchedQueries/*"
    shell: 
        '''
        mkdir output
        python3 !{projectDir}/File_Compatible.py !{inputFile}
        mv output batchedQueries
        '''
}

process Query{
    input:
        path batchFile
    output:
        path "${batchFile}_output.csv"
    shell:
        '''
        
        query_list=$(cat !{batchFile})

        mysql --batch --user=genome --host=genome-mysql.cse.ucsc.edu -N -A -D hg38 -e \
        'select name,
        name2,
        chrom, 
        wg.strand,
        txStart,
        txEnd,
        exonCount,
        exonStarts,
        exonEnds,
        genoStart,
        genoEnd, 
        repName, 
        repFamily, 
        repClass,
        rmsk.strand,
        cdsStart,
        cdsEnd
        from wgEncodeGencodeBasicV39 as wg, rmsk 
        where name2 in'$query_list'and 
        wg.chrom = rmsk.genoName and 
        cast(wg.txStart as signed)-1500 <= cast(rmsk.genoStart as signed)and
        cast(wg.txEnd as signed)+1500 >= cast(rmsk.genoEnd as signed);'\
        > !{batchFile}_sql_out.tsv &&\
        python3 -u !{projectDir}/OverwriTE_V13.py -SQL !{batchFile}_sql_out.tsv -Output !{batchFile}_output.csv 
        
        '''
}

process Gather{
    input:
        tuple file(batchedOutput)
        val outputFile
    output: 
        path "${outputFile}"
    shell: 
        '''
        echo "!{batchedOutput}"
        echo 'name,name2,chrom,genoStrand,genoLength,repName,repStart,repLength,repFamily,repClass,repStrand,classification,len_classification,regionLength,compleTE_seq' > !{outputFile}
        cat !{batchedOutput}|grep -v name,name2,chrom >> !{outputFile}

        '''
}