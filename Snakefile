

output_dir=config['output_dir']
sample_file=output_dir + config['sample_file']
gtf_file=output_dir + 'data/gtfs/all_tissues.combined.gtf'
rule all:
    input: 'results.html'


rule overall_stats:
    input: gtf_file
    output:'data/overall_stats.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/overall_stats.R {output_dir} {sample_file} {input} {output}
        '''

rule novel_exon_catagories:
    input: output_dir + 'rdata/novel_exon_classification.rdata', gtf_file
    output: 'data/exon_classification.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/novel_exon_stat_info.R {output_dir} {input} {output}
        '''
rule number_of_tx_per_tissue:
    input:gtf_file, sample_file, output_dir+'data/seqs/pep_fasta_meta_info.tsv', output_dir+'rdata/novel_exon_classification.rdata',\
    output_dir+'testing/transcripts.fa.transdecoder.genome.gff3', output_dir+'data/misc/gfc_TCONS_to_st_MSTRG.tsv'
    output:'data/novel_tx_by_tissue.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/number_of_tx_per_tissue.R {output_dir} {input} {output}
        '''

rule splicing_hm:
    input:output_dir+'data/rmats/all_tissues_psi.tsv', sample_file
    outpur:'data/splicing_heatmap.Rdata'
    shell:
        '''
        module load R
        Rscript scripts/splicing_heatmap.R {output_dir} {input} {output}
        '''

rule knit_notebooks:
    input:'data/overall_stats.Rdata', 'data/exon_classification.Rdata', 'data/novel_tx_by_tissue.Rdata','data/splicing_heatmap.Rdata'
    output: 'results_v1.html'
    shell:
        '''
        module load R
        Rscript ~/scripts/render_rmd.R notebooks/results_v1.Rmd
        '''
