process Overwrite_wrapper { 
cpus 6
input: 
    file gene_names
    file file_Compatible_py
    file overwritev5
    file overwritev11
output: 
    file ("OverwriTE_annotations.csv")


shell:
'''
#!/bin/bash

mkdir output

python3 !{file_Compatible_py} !{gene_names}

start=1
end=21

for (( c=$start; c<=$end; c++ ))
do
	echo 'Commencing Batch'$c
	## touch output/Batch$c
	./!{overwritev5} output/Batch$c 
done

echo 'name,name2,chrom,genoStrand,genoLength,repName,repStart,repLength,repFamily,repClass,repStrand,classification,len_classification' > output/OverwriTE_annotations.csv

cat output/Batch?_output.csv|grep -v name,name2,chrom >> OverwriTE_annotations.csv


'''
}