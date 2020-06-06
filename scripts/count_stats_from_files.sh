#!/bin/bash
sample_table=$1
data_dir=$2
dntx_quant_path=$3
DNTX_mapping_rates=$4
gencode_mappingrates=$5
gtfs_counts=$6
enst_tx_counts=$7
union_enst_tx=$8
rm -f $4 $5 $6 $7 $8
subtissues=`cut ${sample_table} -f4 | tail -n+2  |sort -u `

for subt in $subtissues
do
    for sample in ${dntx_quant_path}/${subt}/*/aux_info/meta_info.json
        do
            echo $sample `grep num_processed $sample ` `grep percent_mapped $sample` >> ${DNTX_mapping_rates}
        done
done


for sample in ${dntx_quant_path}/gencode/*/aux_info/meta_info.json
    do
        echo $sample `grep num_processed $sample ` `grep percent_mapped $sample` >> ${gencode_mappingrates}
    done


for gtf in ${data_dir}/data/gtfs/*/*.gtf
    do
        echo $gtf `awk '$3 == "transcript"' $gtf | wc -l` >> ${gtfs_counts}
    done

for gtf in ${data_dir}/data/gtfs/raw_tissue_gtfs/*.combined.gtf
    do
        distinct_enst=`grep -e 'class_code \"=\"' $gtf  | grep -o -P 'cmp_ref \"\w+\.\d\"' - | cut -d'"' -f2 | sort -u`
        echo ${distinct_enst}|tr ' ' '\n' >> ${union_enst_tx}
        echo $gtf `echo ${distinct_enst}|tr ' ' '\n' | wc -l` >> ${enst_tx_counts}
    done 