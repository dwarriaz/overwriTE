#!/bin/bash

mkdir output

python3 File_Compatible.py $1

start=1
end=21

for (( c=$start; c<=$end; c++ ))
do
	echo 'Commencing Batch'$c
	## touch output/Batch$c
	./OverwriTE_V5.sh output/Batch$c 
done

echo 'name,name2,chrom,genoStrand,genoLength,repName,repStart,repLength,repFamily,repClass,repStrand,classification,len_classification' > output/OverwriTE_annotations.csv

cat output/Batch*_output.csv|grep -v name,name2,chrom >> output/OverwriTE_annotations.csv

