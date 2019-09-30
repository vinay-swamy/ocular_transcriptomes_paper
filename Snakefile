'''
Snakemake wrapper to generate data for occular_transcriptomes_paper

'''
working_dir=config['working_dir']
output_dir=config['output_dir']
sample_file=output_dir + config['sample_file']
gtf_file=output_dir + 'data/gtfs/all_tissues.combined.gtf'
rule all:
    input: 'notebooks/results_v1.html'

rule transcriptome_pipeline_strategies:
    input:  output_dir + 'data/gtfs/raw_tissue_gtfs/RPE_Fetal.Tissue_st.gtf'
    output: 'clean_data/lib_sizes.tab'
    shell:
        '''
        bash scripts/count_stats_from_files.sh
        module load gffcompare
        gffcompare -r {output_dir}/ref/gencode_comp_ano.gtf -o /tmp/tissue {input}
        module load
        Rscript scripts/transcriptome_pipe_line_stats.R


        '''


rule overall_stats:
    input: gtf_file, output_dir + 'rdata/novel_exon_classification.rdata'
    output:working_dir+'clean_data/overall_stats.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/overall_stats.R {output_dir} {input} {sample_file} {output}
        '''

rule novel_exon_catagories:
    input: output_dir + 'rdata/novel_exon_classification.rdata', gtf_file
    output: working_dir+'clean_data/exon_classification.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/novel_exon_stats_info.R  {output_dir} {input} {output}
        '''
rule number_of_tx_per_tissue:
    input:gtf_file, sample_file, output_dir+'data/seqs/pep_fasta_meta_info.tsv',\
    output_dir+'rdata/novel_exon_classification.rdata',\
    output_dir+'data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3',\
    output_dir+'data/misc/gfc_TCONS_to_st_MSTRG.tsv'
    output:working_dir+'clean_data/novel_tx_by_tissue.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/number_of_tx_per_tissue.R {output_dir} {input} {output}
        '''

rule splicing_hm:
    input:output_dir+'data/rmats/all_tissues_psi.tsv', sample_file
    output:working_dir+'clean_data/splicing_heatmap.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/splicing_heatmap.R {output_dir} {input} {output}
        '''

rule knit_notebooks:
    input:working_dir+'clean_data/overall_stats.Rdata',\
    working_dir+'clean_data/exon_classification.Rdata',\
    working_dir+'clean_data/novel_tx_by_tissue.Rdata',\
    working_dir+'clean_data/splicing_heatmap.Rdata'
    output: 'notebooks/results_v1.html'
    shell:
        '''
        module load R
        Rscript ~/scripts/render_rmd.R notebooks/results_v1.Rmd
        '''
