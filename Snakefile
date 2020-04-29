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
tcons2mstrg=data_dir + 'data/gtfs/all_tissues.convtab'

eye_tissues=['Retina_Adult.Tissue', 'Retina_Fetal.Tissue', 'RPE_Adult.Tissue', 'RPE_Fetal.Tissue', 'Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue']

rule all:
    input: 'notebooks/results_v2.html'

rule transcriptome_pipeline_stats:
    input:  
        t2m=tcons2mstrg,
    params: 
        tissue='RPE_Fetal.Tissue'
    output: 
        DNTX_mr= working_dir + 'clean_data/DNTX_salmon_mapping_rates.tab',\
        gencode_mr= working_dir + 'clean_data/gencode_salmon_mapping_rates.tab',\
        tx_counts= working_dir+'clean_data/all_gtfs_tx_counts.tab',\
        clean_data= working_dir+'clean_data/rdata/transcriptome_pipeline_stats.Rdata'
    shell:
        '''
        bash scripts/count_stats_from_files.sh \
            {data_dir} \
            {output.DNTX_mr} \
            {output.gencode_mr} \
            {output.tx_counts}
        module load {R_version}
        Rscript scripts/transcriptome_pipeline_stats.R \
            --workingDir {working_dir}\
            --dataDir {data_dir}\
            --sampleTableFile {sample_file}\
            --anoGtfFile {gtf_file}\
            --tcons2mstrgFile {input.t2m}\
            --tissue {params.tissue}\
            --gencodeMapRateFile {output.gencode_mr}\
            --DNTXMapRateFile {output.DNTX_mr}\
            --gtfTxCountFile {output.tx_counts}\
            --cleanData {output.clean_data}
        '''

rule summarizeBuildResults:
    input:
        exon_class=exon_classification_file , 
        tc2mstrg=tcons2mstrg ,gff3=gff3_file  
    output:
        color_df=working_dir+ 'clean_data/rdata/tissue_to_colors.Rdata', 
        plotting_data= working_dir+ 'clean_data/rdata/buildResultsSummary.Rdata'
    shell:
        '''
        module load {bedtools_version}
        module load {R_version}
        which R 
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
    input:
        eye_gtfs=expand( data_dir+'data/gtfs/final_tissue_gtfs/{tissue}.gtf', tissue=eye_tissues), 
        pan_body_gtf=gtf_file, 
        DNTX_mr= working_dir + 'clean_data/DNTX_salmon_mapping_rates.tab',
        gencode_mr= working_dir + 'clean_data/gencode_salmon_mapping_rates.tab',
    params: 
        final_gtf_dir = data_dir + 'data/gtfs/final_gtfs/', 
        raw_gtf_dir = data_dir + 'data/gtfs/raw_tissue_gtfs/', 
        ref_tx_exon_rdata = data_dir + 'rdata/all_ref_tx_exons.rdata', 
        core_tight_file = data_dir + 'ref/core_tight.Rdata'
    output: 
        pan_eye_gtf = 'clean_data/pan_eye_txome.combined.gtf', 
        clean_data = working_dir + 'clean_data/rdata/paper_numbers_and_sup_figs.Rdata'
    shell:
        '''
        module load gffcompare 
        module load {R_version}
        gffcompare -r {ref_gtf} --strict-match -o clean_data/pan_eye_txome {input.eye_gtfs} 
        Rscript scripts/paper_numbers_sup_figs.R \
            --workingDir {working_dir} \
            --dataDir {data_dir} \
            --sampleTableFile {sample_file} \
            --allTissueGtf {input.pan_body_gtf} \
            --EyeOnlyGtf {output.pan_eye_gtf} \
            --pathToRawGtfs {params.raw_gtf_dir} \
            --allRefAno {params.ref_tx_exon_rdata} \
            --coreTight {params.core_tight_file} \
            --outNumFile {output.clean_data} \
            --dntxMapRate {input.DNTX_mr} \
            --gencodeMapRate {input.gencode_mr}

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
            --quantFile {input.quant} \
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
    working_dir + 'clean_data/rdata/paper_numbers_and_sup_figs.Rdata' 
    #working_dir + 'clean_data/rdata/fetal_novel_de_results.Rdata', \
    #working_dir + 'clean_data/rdata/fetal_novel_de_hm.Rdata',\
    output: 'notebooks/results_v2.html'
    shell:
        '''
        module load {R_version}
        Rscript ~/scripts/render_rmd.R notebooks/results_v2.Rmd
        '''

