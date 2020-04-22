#!/bin/bash
quant_path=$1
DNTX_mapping_rates=$2
gencode_mappingrates=$3
gtfs_counts=$4
rm -f $2 $3 $4
mkdir -p tmp
for sample in ${quant_path}/data/salmon_quant/dntx/*/*/aux_info/meta_info.json
    do
        echo $sample `grep num_processed $sample ` `grep percent_mapped $sample` >> ${DNTX_mapping_rates}
    done


for sample in ${quant_path}/gencode_quant/quant_files/*/aux_info/meta_info.json
    do
        echo $sample `grep num_processed $sample ` `grep percent_mapped $sample` >> ${gencode_mappingrates}
    done


for gtf in ${quant_path}/data/gtfs/*/*.gtf
    do
        echo $gtf `awk '$3 == "transcript"' $gtf | wc -l` >> ${gtfs_counts}
    done
