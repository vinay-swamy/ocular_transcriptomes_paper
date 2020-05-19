'''
Snakemake wrapper to generate data for occular_transcriptomes_paper

'''
import yaml
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
with open(config['subtissue_file']) as stfl:
    subtissues = [line.strip('\n') for line in stfl ]
working_dir=config['working_dir']
longread_dir = config['longread_dir']
data_dir=config['data_dir']
fql=config['fastq_path']
sample_file=config['sample_file']
files_yaml = config['files_yaml']
with open(config['files_yaml']) as fyml:
    files=yaml.load(fyml,Loader=yaml.FullLoader)
ref_gtf = files['ref_gtf']

eye_tissues=['Retina_Adult.Tissue', 'Retina_Fetal.Tissue', 'RPE_Adult.Tissue', 'RPE_Fetal.Tissue', 'Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue']

rule all:
    input: 'notebooks/results_v2.html'

rule calc_txome_tx_counts_and_mapping_Rates:
    input:  
        full_anno_gtf =files['anno_gtf'],
        tcons2mstrg = files['tcons2mstrg']
    output: 
        DNTX_mr= files['DNTX_mr'],
        gencode_mr= files['gencode_mr'],
        gtf_tx_counts= files['gtf_tx_counts'],
        txcounts_mr_rdata= files['txome_stats_rdata']
    shell:
        '''
        bash scripts/count_stats_from_files.sh \
            {sample_file} \
            {data_dir} \
            {output.DNTX_mr} \
            {output.gencode_mr} \
            {output.gtf_tx_counts}
             
        module load {R_version}
        Rscript scripts/transcriptome_pipeline_stats.R \
            --workingDir {working_dir}\
            --dataDir {data_dir}\
            --filesYaml {files_yaml}
        '''

rule summarizeBuildResults:
    input:
        exon_class=files['exon_class_rdata'],
        tc2mstrg=files['tcons2mstrg'],
        full_anno_gtf =files['anno_gtf']
    params: 
        cds_gtf_dir = 'clean_data/CDS_gtf/', 
        raw_cds_track = '/data/swamyvs/ocular_transcriptomes_paper/clean_data/CDS_gtf/CDS_comp.tracking'
    output:
        color_mapping_rdata=files['color_mapping_rdata'],
        build_results_rdata= files['build_results_rdata']
    shell:
        '''
        module load gffcompare
        rm -rf {params.cds_gtf_dir}
        mkdir -p {params.cds_gtf_dir}
        awk '$3 == "CDS"' {ref_gtf} > {params.cds_gtf_dir}/gencode_CDS.gtf
        awk '$3 == "CDS"' {input.full_anno_gtf} > {params.cds_gtf_dir}/alltissuegtf_cds_only.gtf
        gffcompare --strict-match -p CDSID -o {params.cds_gtf_dir}/CDS_distinct {params.cds_gtf_dir}/alltissuegtf_cds_only.gtf {params.cds_gtf_dir}/gencode_CDS.gtf
        gffcompare --strict-match -r {params.cds_gtf_dir}CDS_distinct.combined.gtf -o {params.cds_gtf_dir}/CDS_comp {params.cds_gtf_dir}/alltissuegtf_cds_only.gtf   

        module load {bedtools_version}
        module load {R_version}
        Rscript scripts/summarize_build_results.R \
            --workingDir {working_dir}\
            --dataDir {data_dir}\
            --rawCDStrack {params.raw_cds_track} \
            --filesYaml {files_yaml}
   
        '''

rule summarize_long_read_results:
    output: 
        working_dir + 'clean_data/rdata/longread_summary.Rdata'
    shell:
        '''
        module load {R_version}
        Rscript scripts/summarise_longread_results.R \
            --longReadDir {longread_dir} \
            --cleanedData {output}
        '''

rule paper_numbers_and_sup_figs:
    input:
        pan_eye_gtf = files['pan_eye_gtf'],
        all_tissue_gtf=files['anno_gtf'],
        DNTX_mr= files['DNTX_mr'],
        gencode_mr= files['gencode_mr']
    output: 
        clean_data = files['paper_numbers_rdata']
    shell:
        '''
        module load {R_version}
        Rscript scripts/paper_numbers_sup_figs.R \
            --workingDir {working_dir} \
            --dataDir {data_dir} \
            --filesYaml {files_yaml}
        '''


rule novel_isoforms_ocular_tissues:
    input: 
        quant = files['all_tissue_quant'], 
        exon_class_rdata = files['exon_class_rdata'],
        gtf_ano_file = files['anno_gtf']
    output: 
        outdata= files['novel_isoform_analysis_rdata']
    shell:
        '''
        module load {R_version}
        Rscript scripts/novel_isoforms_ocular_tissues.R \
            --workingDir {working_dir} \
            --dataDir {data_dir} \
            --filesYaml {files_yaml}
            
        '''

rule process_VEP:
    # input: 
    #     expand(data_dir + 'data/vep/{subtissue}/variant_summary.txt', subtissue = subtissues)
    output: 
        all_variant_results = working_dir + 'clean_data/rdata/vep_all_alleles.Rdata',
        example_variant_results = working_dir +'clean_data/rdata/vep_eye_example.Rdata'
    shell:
        '''
        module load {R_version}
        Rscript scripts/process_VEP.R \
        --dataDir {data_dir} \
        --allVariantFile   {output.all_variant_results} \
        --exampleVariantFile {output.example_variant_results}     
        '''


rule process_hmmer:
    input:
        data_dir + 'data/novel_loci/hmmer/seq_hits.tsv'
    output:
        working_dir + 'clean_data/rdata/hmmer_results.Rdata'
    shell:
        '''
        module load {R_version}
        Rscript scripts/process_HMMER.R \
            --dataDir {data_dir} \
            --filesYaml {files_yaml}
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
    input: 
        working_dir + 'clean_data/rdata/buildResultsSummary.Rdata', 
        working_dir + 'clean_data/rdata/transcriptome_pipeline_stats.Rdata', 
        working_dir + 'clean_data/rdata/novel_isoforms.Rdata', 
        working_dir + 'clean_data/rdata/paper_numbers_and_sup_figs.Rdata', 
        working_dir + 'clean_data/rdata/longread_summary.Rdata',
        working_dir +'clean_data/rdata/vep_eye_example.Rdata', 
        working_dir + 'clean_data/rdata/hmmer_results.Rdata'
    #working_dir + 'clean_data/rdata/fetal_novel_de_results.Rdata', \
    #working_dir + 'clean_data/rdata/fetal_novel_de_hm.Rdata',\
    output: 
        'notebooks/results_v2.html'
    shell:
        '''
        module load {R_version}
        Rscript ~/scripts/render_rmd.R notebooks/results_v2.Rmd
        '''

