'''
Snakemake wrapper to generate data for occular_transcriptomes_paper

'''

def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        file.readline()#skip the first line because it has a header
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'paired':True if info[1]=='y' else False, 'tissue':info[2],'subtissue':info[3]}
    return(res)

def salmon_input(id,sample_dict,fql):
    paired=sample_dict[id]['paired']
    id= fql + 'fastq_files/' + id
    if paired:
        return('-1 {s}_1.fastq.gz -2 {s}_2.fastq.gz'.format(s=id))
    else:
        return('-r {}.fastq.gz'.format(id))

def subtissue_to_sample(subtissue, sample_dict):
    res=[]
    [res.append(f'salmon_experiment/quant/{subtissue}/{sample}/quant.sf') for sample in sample_dict.keys() if sample_dict[sample]['subtissue'] == subtissue ]
    return(res)


salmon_version= config['salmon_version']
R_version= config['R_version']
bedtools_version= config['bedtools_version'] 
TransDecoder_version=config['TransDecoder_version']

# sd_salmonexp= readSampleFile('salmon_experiment/sampletabSE.tsv')

working_dir=config['working_dir']
data_dir=config['data_dir']
fql=config['fastq_path']
sample_file=config['sample_file']
gtf_file=data_dir + 'data/gtfs/all_tissues.combined.gtf'
exon_classification_file=data_dir + 'data/rdata/novel_exon_classification.Rdata'
gff3_file=data_dir + 'data/seqs/transdecoder_results/all_tissues.combined_transdecoderCDS.gff3'
ref_gtf= data_dir + 'ref/gencode_comp_ano.gtf'
tcons2mstrg=data_dir + 'data/misc/TCONS2MSTRG.tsv'

eye_tissues=['Retina_Adult.Tissue', 'Retina_Fetal.Tissue', 'RPE_Adult.Tissue', 'RPE_Fetal.Tissue', 'Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue']

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
        module load {R_version}
        Rscript scripts/transcriptome_pipeline_stats.R \
            --workingDir {working_dir}\
            --dataDir {data_dir}\
            --sampleTableFile {sample_file}\
            --anoGtfFile {gtf_file}\
            --pathToBs {params.path_to_bs}\
            --tcons2mstrgFile {input.t2m}\
            --tissue {params.tissue}\
            --gencodeMapRateFile {output.gencode_mr}\
            --DNTXMapRateFile {output.DNTX_mr}\
            --gtfTxCountFile {output.tx_counts}\
            --cleanData {output.clean_data}
        '''

rule summarizeBuildResults:
    input:exon_class=exon_classification_file , tc2mstrg=tcons2mstrg ,gff3=gff3_file  
    output:color_df=working_dir+ 'clean_data/rdata/tissue_to_colors.Rdata', plotting_data= working_dir+ 'clean_data/rdata/buildResultsSummary.Rdata'
    shell:
        '''
        module load {bedtools_version}
        module load {R_version}
        Rscript scripts/summarize_build_results.R \
            --workingDir {working_dir}\
            --dataDir {data_dir}\
            --exonClassFile {input.exon_class}\
            --gtfFile {gtf_file}\
            --sampleTableFile {sample_file}\
            --gff3File {input.gff3}\
            --tcons2mstrgFile {input.tc2mstrg}\
            --colorMappingDf {output.color_df}\
            --dataToPlot {output.plotting_data}
        '''


rule paper_numbers_and_sup_figs:
    input:eye_gtfs=expand( data_dir+'data/gtfs/final_gtfs/{tissue}.gtf', tissue=eye_tissues), pan_body_gtf=gtf_file 
    params: final_gtf_dir = data_dir + 'data/gtfs/final_gtfs/', raw_gtf_dir = data_dir + 'data/gtfs/raw_tissue_gtfs/', ref_tx_exon_rdata = data_dir + 'rdata/all_ref_tx_exons.rdata', core_tight_file = data_dir + 'ref/core_tight.Rdata',
    output: pan_eye_gtf = 'clean_data/pan_eye_txome.combined.gtf', clean_data = working_dir + 'clean_data/rdata/paper_numbers_and_sup_figs.Rdata'
    shell:
        '''
        module load gffcompare 
        module load {R_version}
        gffcompare -r {ref_gtf} --strict-match -o clean_data/pan_eye_txome {input.eye_gtfs} 
        Rscript scripts/paper_numbers_sup_figs.R \
            {working_dir} \
            {data_dir} \
            {sample_file} \
            {input.pan_body_gtf} \
            {output.pan_eye_gtf} \
            {params.raw_gtf_dir} \
            {params.ref_tx_exon_rdata} \
            {params.core_tight_file} \
            {output.clean_data}

        '''

rule splicing_hm:
    input:psi_file=data_dir+'data/rmats/all_tissues_psi.tsv', tc2mstrg=tcons2mstrg, classfile=exon_classification_file, color_df=working_dir + 'clean_data/rdata/tissue_to_colors.Rdata', gtf_ano_file = data_dir + 'data/gtfs/all_tissues.combined_NovelAno.gtf'
    output:working_dir+ 'clean_data/plots/novel_exon_splicing_heatmap.tiff'
    shell:
        '''
        module load {R_version}
        Rscript scripts/splicing_heatmap.R \
            --workingDir {working_dir} \
            --psiFile {input.psi_file} \
            --sampleFile {sample_file} \
            --tc2mstrgFile {input.tc2mstrg} \
            --exonClassFile {input.classfile} \
            --gtfFile {input.gtf_ano_file} \
            --colorDf {input.color_df} \
            --outFile {output}
        '''


rule novel_isoforms_ocular_tissues:
    input: quant = data_dir + 'data/all_tissue_quant.Rdata', gtf_ano_file = data_dir + 'data/gtfs/all_tissues.combined_NovelAno.gtf'
    output: outdata= working_dir + 'clean_data/rdata/novel_isoforms.Rdata'
    shell:
        '''
        module load {R_version}
        Rscript scripts/novel_isoforms_ocular_tissues.R \
            --workingDir {working_dir} \
            --dataDir {data_dir} \
            --novelExonClassFile {exon_classification_file} \
            --gtfFile {input.gtf_ano_file} \
            --sampleTableFile {sample_file} \
            --tcons2mstrgFile {tcons2mstrg} \
            --quanFile {input.quant} \
            --cleanData {output.outdata}
            
        '''

# rule novel_tx_in_fetal_retina_analysis:
#     input:eiad='clean_data/EiaD_quant.Rdata', tc2mstrg=tcons2mstrg , gff3=gff3_file
#     output:working_dir+ 'clean_data/rdata/fetal_novel_de_results.Rdata',working_dir+ 'clean_data/fetal_de_novel_tx.txt', working_dir+ 'clean_data/rdata/fetal_novel_de_hm.Rdata'
#     shell:
#         '''
#         module load {R_version}
#         Rscript scripts/analyze_novel_tx_fetal_retina.R {working_dir} {sample_file} {input.eiad} {input.tc2mstrg} {gtf_file} {input.gff3} {output}
#         '''

'''
Get Nv tx >  remove from fasta > build index
run this to start
tail -n+2  /data/swamyvs/eyeintegration_splicing/data/salmon_quant/Retina_Fetal.Tissue/SRS897012/quant.sf | grep Retina_Fetal.tissue - | cut -f1 > clean_data/salmon_experiment/Retina_Fetal.Tissue_tx_to_remove.tab

'''


rule knit_notebooks:
    input: \
    working_dir + 'clean_data/rdata/buildResultsSummary.Rdata', \
    working_dir + 'clean_data/rdata/transcriptome_pipeline_stats.Rdata', \
    working_dir + 'clean_data/rdata/novel_isoforms.Rdata', \
    working_dir + 'clean_data/plots/novel_exon_splicing_heatmap.tiff', \
    working_dir + 'clean_data/rdata/paper_numbers_and_sup_figs.Rdata' 
    #working_dir + 'clean_data/rdata/fetal_novel_de_results.Rdata', \
    #working_dir + 'clean_data/rdata/fetal_novel_de_hm.Rdata',\
    output: 'notebooks/results_v2.html'
    shell:
        '''
        module load {R_version}
        Rscript ~/scripts/render_rmd.R notebooks/results_v2.Rmd
        '''



# rule setup_salmon_exp:
#     output: lens='salmon_experiment/Lens_Stem.Cell.Line_tx_to_remove.tab', ret='salmon_experiment/Retina_Fetal.Tissue_tx_to_remove.tab', st= 'salmon_experiment/sampletabSE.tsv'
#     shell:
#         '''
#         tail -n+2  /data/swamyvs/eyeintegration_splicing/data/salmon_quant/Retina_Fetal.Tissue/SRS897012/quant.sf | grep Retina_Fetal.Tissue - | cut -f1 > salmon_experiment/Retina_Fetal.Tissue_tx_to_remove.tab

#         tail -n+2   /data/swamyvs/eyeintegration_splicing/data/salmon_quant/Lens_Stem.Cell.Line/SRS1747747/quant.sf | grep  Lens_Stem.Cell.Line  - | cut -f1 > salmon_experiment/Lens_Stem.Cell.Line_tx_to_remove.tab

#         grep Retina_Fetal.Tissue sampleTableFull.tsv > salmon_experiment/sampletabSE.tsv
#         grep Lens_Stem.Cell.Line sampleTableFull.tsv >> salmon_experiment/sampletabSE.tsv

#         '''


# rule build_salmon_index_no_novel:
#     input:tx_to_remove='salmon_experiment/{tissue}_tx_to_remove.tab', fasta= data_dir + 'data/seqs/{tissue}_tx.fa'
#     output:outfasta='salmon_experiment/{tissue}.fa', index= directory('salmon_experiment/{tissue}')
#     shell:
#         ''' 
#         python3 scripts/remove_entries_from_fasta.py {input.fasta} {input.tx_to_remove} {output.outfasta}
#         module load {salmon_version}
#         salmon index -t {output.outfasta} -i {output.index} --type quasi --perfectHash -k 31
#         '''

# rule run_salmon:
#     input: fastqs=lambda wildcards: [fql+f'fastq_files/{wildcards.sampleID}_1.fastq.gz',fql+f'fastq_files/{wildcards.sampleID}_2.fastq.gz'] if sd_salmonexp[wildcards.sampleID]['paired'] else fql+f'fastq_files/{wildcards.sampleID}.fastq.gz',\
#         index='salmon_experiment/{tissue}'
#     params: cmd=lambda wildcards: salmon_input(wildcards.sampleID,sd_salmonexp,fql),\
#      outdir=lambda wildcards: f'salmon_experiment/quant/{wildcards.tissue}/{wildcards.sampleID}'
#     output: 'salmon_experiment/quant/{tissue}/{sampleID}/quant.sf'
#     shell:
#         '''
#         id={wildcards.sampleID}
#         module load {salmon_version}
#         salmon quant -p 8 -i {input.index} -l A --gcBias --seqBias --validateMappings {params.cmd} -o {params.outdir}
#         '''



# rule ORF_exp:
#     input:gtf= data_dir + 'data/gtfs/raw_tissue_gtfs/Retina_Fetal.Tissue.combined.gtf'
#     output:'trandecode_exp/tr_dec/all_Retina_Fetal.Tissue_tx.gff3'
#     shell:
#         '''
#         rm -rf TransDecoder
#         git clone https://github.com/TransDecoder/TransDecoder.git
#         cd TransDecoder
#         module load {TransDecoder_version}
#         mkdir -p ../trandecode_exp/tr_dec
#         ./util/gtf_genome_to_cdna_fasta.pl {input.gtf} ../ref/gencode_genome.fa > transcripts.fasta
#         ./util/gtf_to_alignment_gff3.pl {input.gtf}> transcripts.gff3
#         TransDecoder.LongOrfs -t tran scripts.fasta
#         TransDecoder.Predict --single_best_only -t transcripts.fasta
#         ./util/cdna_alignment_orf_to_genome_orf.pl \
#             transcripts.fasta.transdecoder.gff3 \
#             transcripts.gff3 \
#             transcripts.fasta > ../trandecode_exp/tr_dec/all_Retina_Fetal.Tissue_tx.gff3
#         '''


# rule gffcompare_all:
#     output:'gffcomp_test/all_samples.combined.gtf'
#     shell:
#         '''
#         module load gffcompare 
#         input=`ls /data/swamyvs/eyeintegration_splicing/st_out/*.gtf | tr '\\n' ' '`
#         gffcompare -r /data/swamyvs/ref/gencode_comp_ano.gtf -o gffcomp_test/all_samples $input

#         '''

