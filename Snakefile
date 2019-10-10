'''
Snakemake wrapper to generate data for occular_transcriptomes_paper

'''
working_dir=config['working_dir']
output_dir=config['output_dir']
sample_file=config['sample_file']
gtf_file=output_dir + 'data/gtfs/all_tissues.combined.gtf'
exon_classification_file=output_dir + 'data/rdata/novel_exon_classification.Rdata'
gff3_file=output_dir + 'data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3'
tcons2mstrg=output_dir + 'data/misc/TCONS2MSTRG.tsv'

rule all:
    input: 'notebooks/results_v1.html'

rule transcriptome_pipeline_stats:
    input:  t2m=tcons2mstrg,
    params: path_to_bs=output_dir+'data/rawST_tx_quant_files/RPE_Fetal.Tissue/'
    output: 'clean_data/transcriptome_pipeline_stats.Rdata'
    shell:
        '''
        bash scripts/count_stats_from_files.sh {output_dir}
        module load R
        Rscript scripts/transcriptome_pipeline_stats.R {working_dir} {sample_file} {input.t2m} {gtf_file} {params.path_to_bs} {output}
        '''

rule summarizeBuildResults:
    input:exon_class=exon_classification_file , gff3=gff3_file, tc2mstrg=tcons2mstrg
    output:color_df='clean_data/tissue_to_colors.Rdata', plotting_data='clean_data/buildResultsSummary.Rdata'
    shell:
        '''
        module load bedtools
        module load R
        Rscript scripts/summarize_build_results.R {working_dir} {input.exon_class} {gtf_file} {sample_file} {input.gff3} {input.tc2mstrg} {output.color_df} {output.plotting_data}
        '''

# rule overall_stats:
#     input: gtf_file, output_dir + 'rdata/novel_exon_classification.rdata'
#     output:working_dir+'clean_data/overall_stats.Rdata'
#     shell:
#         '''
#         module load R
#         Rscript scripts/overall_stats.R {output_dir} {input} {sample_file} {output}
#         '''
#
# rule novel_exon_catagories:
#     input: output_dir + 'rdata/novel_exon_classification.rdata', gtf_file
#     output: working_dir+'clean_data/exon_classification.Rdata'
#     shell:
#         '''
#         module load R
#         Rscript scripts/novel_exon_stats_info.R  {output_dir} {input} {output}
#         '''
# rule number_of_tx_per_tissue:
#     input:gtf_file, sample_file, output_dir+'data/seqs/pep_fasta_meta_info.tsv',\
#     output_dir+'rdata/novel_exon_classification.rdata',\
#     output_dir+'data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3',\
#     output_dir+'data/misc/gfc_TCONS_to_st_MSTRG.tsv'
#     output:working_dir+'clean_data/novel_tx_by_tissue.Rdata'
#     shell:
#         '''
#         module load R
#         Rscript scripts/number_of_tx_per_tissue.R {output_dir} {input} {output}
#        '''

rule splicing_hm:
    input:psi_file=output_dir+'data/rmats/all_tissues_psi.tsv', tc2mstrg=tcons2mstrg, classfile=exon_classification_file, color_df='clean_data/tissue_to_colors.Rdata'
    output:'clean_data/splicing_analysis_results.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/splicing_heatmap.R {working_dir} {input.psi_file} {sample_file} {input.tc2mstrg} {input.classfile} {gtf_file} {input.color_df} {output}
        '''

rule novel_tx_in_fetal_retina_analysis:
    input:eiad='clean_data/EiaD_quant.Rdata', tc2mstrg=tcons2mstrg, gff3=gff3_file
    output:'clean_data/fetal_novel_de_results.Rdata','clean_data/fetal_de_novel_tx.txt', 'clean_data/fetal_novel_de_hm.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/analyze_novel_tx_fetal_retina.R {working_dir} {sample_file} {input.eiad} {input.tc2mstrg} {gtf_file} {input.gff3} {output}
        '''



rule knit_notebooks:
    input:'clean_data/buildResultsSummary.Rdata', 'clean_data/splicing_analysis_results.Rdata', 'clean_data/fetal_novel_de_results.Rdata', 'clean_data/fetal_novel_de_hm.Rdata','clean_data/transcriptome_build_summary.Rdata'

    output: 'notebooks/results_v1.html'
    shell:
        '''
        module load R
        Rscript ~/scripts/render_rmd.R notebooks/results_v1.Rmd
        '''
