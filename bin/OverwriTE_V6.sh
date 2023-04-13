#!/bin/bash 

query_list=$(cat $1)

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
    cdsEnds
    from wgEncodeGencodeBasicV41 as wg, rmsk 
    where name2 in'$query_list'and 
	wg.chrom = rmsk.genoName and 
	cast(wg.txStart as signed)-1500 <= cast(rmsk.genoStart as signed)and
	cast(wg.txEnd as signed)+1500 >= cast(rmsk.genoEnd as signed);'\
    > $1_sql_out.tsv &&\
     python3 -u OverwriTE_V11.py -SQL $1_sql_out.tsv -Output $1_output.csv 
