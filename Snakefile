'''
Snakemake wrapper to generate data for occular_transcriptomes_paper

'''
working_dir=config['working_dir']
data_dir=config['data_dir']
sample_file=config['sample_file']
gtf_file=data_dir + 'data/gtfs/all_tissues.combined.gtf'
exon_classification_file=data_dir + 'data/rdata/novel_exon_classification.Rdata'
#gff3_file=data_dir + 'data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3'
gff3_file='tmp/dummy.gff3'
tcons2mstrg=data_dir + 'data/misc/TCONS2MSTRG.tsv'

rule all:
    input: 'notebooks/results_v2.html'

rule transcriptome_pipeline_stats:
    input:  t2m=tcons2mstrg,
    params: path_to_bs=data_dir+'data/salmon_quant/RPE_Fetal.Tissue/', tissue='RPE_Fetal.Tissue'
    output: DNTX_mr= working_dir + 'clean_data/DNTX_salmon_mapping_rates.tab',\
            gencode_mr= working_dir + 'clean_data/gencode_salmon_mapping_rates.tab',\
            tx_counts= working_dir+'clean_data/all_gtfs_tx_counts.tab',\
            clean_data= working_dir+'clean_data/rdata/transcriptome_pipeline_stats.Rdata'
    shell:
        '''
        bash scripts/count_stats_from_files.sh {data_dir} {output.DNTX_mr} {output.gencode_mr} {output.tx_counts}
        module load R
        Rscript scripts/transcriptome_pipeline_stats.R \
            {working_dir}\
            {data_dir}\
            {sample_file}\
            {gtf_file}\
            {params.path_to_bs}\
            {input.t2m}\
            {params.tissue}\
            {output.DNTX_mr}\
            {output.gencode_mr}\
            {output.tx_counts}\
            {output.clean_data}
        '''

rule summarizeBuildResults:
    input:exon_class=exon_classification_file , tc2mstrg=tcons2mstrg ,gff3=gff3_file  
    output:color_df=working_dir+ 'clean_data/rdata/tissue_to_colors.Rdata', plotting_data= working_dir+ 'clean_data/rdata/buildResultsSummary.Rdata'
    shell:
        '''
        module load bedtools
        module load R
        Rscript scripts/summarize_build_results.R \
            {working_dir}\
            {data_dir}\
            {input.exon_class}\
            {gtf_file}\
            {sample_file}\
            {input.gff3}\
            {input.tc2mstrg}\
            {output.color_df}\
            {output.plotting_data}
        '''


rule splicing_hm:
    input:psi_file=data_dir+'data/rmats/all_tissues_psi.tsv', tc2mstrg=tcons2mstrg, classfile=exon_classification_file, color_df=working_dir+ 'clean_data/rdata/tissue_to_colors.Rdata'
    output:working_dir+ 'clean_data/rdata/splicing_analysis_results.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/splicing_heatmap.R {working_dir} {input.psi_file} {sample_file} {input.tc2mstrg} {input.classfile} {gtf_file} {input.color_df} {output}
        '''

rule novel_tx_in_fetal_retina_analysis:
    input:eiad='clean_data/EiaD_quant.Rdata', tc2mstrg=tcons2mstrg , gff3=gff3_file
    output:working_dir+ 'clean_data/rdata/fetal_novel_de_results.Rdata',working_dir+ 'clean_data/fetal_de_novel_tx.txt', working_dir+ 'clean_data/rdata/fetal_novel_de_hm.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/analyze_novel_tx_fetal_retina.R {working_dir} {sample_file} {input.eiad} {input.tc2mstrg} {gtf_file} {input.gff3} {output}
        '''



rule knit_notebooks:
    input:expand(f'{working_dir}{{files}}', files=['clean_data/rdata/buildResultsSummary.Rdata', 'clean_data/rdata/splicing_analysis_results.Rdata', 'clean_data/rdata/fetal_novel_de_results.Rdata', 'clean_data/rdata/fetal_novel_de_hm.Rdata','clean_data/rdata/transcriptome_pipeline_stats.Rdata'])

    output: 'notebooks/results_v2.html'
    shell:
        '''
        module load R
        Rscript ~/scripts/render_rmd.R notebooks/results_v2.Rmd
        '''
