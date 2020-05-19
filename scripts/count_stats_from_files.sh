#!/bin/bash
sample_table=$1
quant_path=$2
DNTX_mapping_rates=$3
gencode_mappingrates=$4
gtfs_counts=$5
rm -f $3 $4 $5
subtissues=`cut ${sample_table} -f4 | sort -u `

for subt in $subtissues
do
    for sample in ${quant_path}/data/salmon_quant/dntx/${subt}/*/aux_info/meta_info.json
        do
            echo $sample `grep num_processed $sample ` `grep percent_mapped $sample` >> ${DNTX_mapping_rates}
        done
done


for sample in ${quant_path}/gencode_quant/quant_files/*/aux_info/meta_info.json
    do
        echo $sample `grep num_processed $sample ` `grep percent_mapped $sample` >> ${gencode_mappingrates}
    done


for gtf in ${quant_path}/data/gtfs/*/*.gtf
    do
        echo $gtf `awk '$3 == "transcript"' $gtf | wc -l` >> ${gtfs_counts}
    done
