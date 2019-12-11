#!/bin/bash
quant_path=$1
rm -rf tmp
mkdir -p tmp
for sample in ${quant_path}/data/salmon_quant/*/*/aux_info/meta_info.json
    do
        echo $sample `grep num_processed $sample ` `grep percent_mapped $sample` >> tmp/DNTX_mappingrates.tab
    done

for sample in ${quant_path}/gencode_quant/quant_files/*/aux_info/meta_info.json
    do
        echo $sample `grep num_processed $sample ` `grep percent_mapped $sample` >> tmp/gencode_mappingrates.tab
    done

for gtf in ${quant_path}/data/gtfs/*/*.gtf
    do
        echo $gtf `awk '$3 == "transcript"' $gtf | wc -l` >> tmp/gtf_tx_counts.tab
    done
