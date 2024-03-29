#!/bin/bash
mkdir output

counter=$(python3 File_Compatible.py $1)
echo $counter
start=1
end=$(echo $counter)

for (( c=$start; c<=$end; c++ ))
do
	echo 'Commencing Batch'$c
	touch output/Batch$c
	./OverwriTE_V5.sh output/Batch$c 
done

echo 'name,name2,chrom,genoStrand,genoLength,repName,repStart,repLength,repFamily,repClass,repStrand,classification,len_classification,regionLength,compleTE_seq' > output/OverwriTE_annotations.csv

cat output/Batch*_output.csv|grep -v name,name2,chrom >> output/OverwriTE_annotations.csv

