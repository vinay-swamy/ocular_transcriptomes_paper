---
title: '*De novo* transcriptomes built from hundreds of eye tissues reveal hundreds of novel gene isoforms'
author: Vinay S Swamy, Temesgen D Fufa, Robert B Hufnagel, David M McGaughey
bibliography: ocular_txome_citations.bib
csl: nucleic-acids-research.csl
output:
  redoc::redoc:
    reference_docx: reference_doc_v1.docx
    keep_md: yes
  html_document:
    df_print: paged
  pdf_document:
    df_print: kable
---

<div class="redoc" id="redoc-codechunk-1">


</div>

<div class="redoc" id="redoc-codechunk-2">


</div>

## Introduction

       The transcriptome is defined as the set of distinct RNA transcripts expressed in a population of identical cells. During transcription several RNA processing steps modify immature RNA and drive the formation of multiple, distinct isoforms for most genes. For example, the human Gencode release 28 contains 97,713 protein coding transcripts across 20,306 genes <span class="redoc" id="redoc-citation-1"><span class="redoc" id="redoc-citation-22"><span class="redoc" id="redoc-citation-25">[-@frankish_gencode_2019]</span></span></span> . RNA processing broadly describes a variety of biological mechanisms and includes alternative promoter usage, alternative splicing, RNA editing, and alternative polyadenylation. The full biological impact of gene isoforms has not been fully elucidated, but multiple studies have shown that gene isoforms can have distinct and critical functions in biological processes like development<span class="redoc" id="redoc-citation-2">[-@dykes_hic2_2018]</span>, cell differentiation<span class="redoc" id="redoc-citation-3"><span class="redoc" id="redoc-citation-11">[-@trapnell_transcript_2010]</span></span>, and cell migration<span class="redoc" id="redoc-citation-4">[-@mitra_splicing_2020]</span>. Alternative usage of isoforms has also been implicated in multiple diseases including cancer<span class="redoc" id="redoc-citation-5">[-@vitting-seerup_landscape_2017]</span>, cardiovascular disease <span class="redoc" id="redoc-citation-6">[-@neagoe_ciprian_titin_2002]</span>, Alzheimer's Disease<span class="redoc" id="redoc-citation-7">[-@mills_rna-seq_2013]</span> and diabetic retinopathy<span class="redoc" id="redoc-citation-8">[-@perrin_diabetic_2005]</span>.(didnt read these papers so need to double check them)

       Some of the first methods using RNA-seq to find novel gene isoforms focused on identifying novel exon-exon junctions and novel exon boundaries based on RNA-seq coverage. <span class="redoc" id="redoc-citation-9">[-@nagalakshmi_transcriptional_2008]</span> More recently, several groups have developed specialized tools to use RNA-seq to reconstruct the whole transcriptome of a biological samples, dubbed *de novo* transcriptome construction <span class="redoc" id="redoc-citation-10"><span class="redoc" id="redoc-citation-36">[-@haas_novo_2013]</span></span>,[-@trapnell_transcript_2010],<span class="redoc" id="redoc-citation-12"><span class="redoc" id="redoc-citation-26">[-@pertea_stringtie_2015]</span></span>.

       *De novo* transcriptome construction uses short (\~100bp) RNA-seq reads to reconstruct full-length mRNA transcripts. However, a large number of samples are necessary to overcome the noise and short read lengths of this type of data. Because of increasingly inexpensive sequencing, data sets of the necessary size are now available. For example, one of the most comprehensive *de novo* transcriptome projects to date has been CHESS, which used the GTEx data set to construct *de novo* transcriptomes in over 9000 RNA-seq samples from 49 distinct location of the body to create a comprehensive annotation of mRNA transcripts across the human body. <span class="redoc" id="redoc-citation-13">[-@gtex_consortium_genetic_2017]</span>,<span class="redoc" id="redoc-citation-14">[-@pertea_chess_2018]</span> However, as the GTEx data set lacks any ocular tissues, the CHESS database is an incomplete annotation of the human transcriptome.

       Despite the increasing number of tools for transcriptome construction there has been no gold standard with which to evaluate precision and sensitivity of *de novo* transcriptome construction on real (not simulated) biological data. Long read sequencing technologies provide a potential solution to this problem as long read sequencing can capture full length transcripts and thus can be used to identify a fuller range of gene isoforms. While long reads have historically been considered inaccurate, the new PacBio Sequel II system sequences long reads as accurately as short read based sequencing <span class="redoc" id="redoc-citation-15">[-@wenger_accurate_2019]</span>.

       We propose that long read based transcriptomes can serve as a ground truth for evaluating short-read base transcriptomes, and in this study we use PacBio long read RNA sequencing to inform the construction of short read transcriptomes. We generated PacBio long read RNA-seq data from a stem cell derived retinal pigmented epithelium (RPE) cell line along with matched Illumina short read RNA-seq. Using the two sources of RNA-seq data we design a rigorous Stringtie-based *de novo* transcriptome pipeline that maximizes the concordance between short and long read *de novo* transcriptomes.

       We apply this pipeline using a previously published data set containing <span class="redoc" id="redoc-citation-16"><span class="redoc" id="redoc-inlinecode-1">368</span> ocular tissue samples compiled from mining publicly available sequencing data <span class="redoc" id="redoc-citation-41"><span class="redoc" id="redoc-citation-68">[-@swamy_eye_2019]</span></span></span>. We use this pipeline to build transcriptomes in three major ocular subtissues: The cornea, retina, and the RPE, using RNA-seq data from both adult and fetal tissues to create a high-quality pan-eye transcriptome. In addition to our ocular samples, we used a subset of the GTEx data set to construct transcriptomes for 49 other locations across the body to facilitate comparisons between transcriptomes across the body.

       We use our gold-standard informed pan eye *de novo* transcriptome to reveal hundreds of novel gene isoforms as well as several novel [ Want to avoid using the word gene ]{.comment-start id="1" author="swamyvs" date="2020-06-23T13:48:29Z"} transcribed genomic loci  []{.comment-end id="1"}in the eye and analyze their potential impact on ocular biology and disease. We provide our *de novo* transcriptomes as a resource to other researchers through an R package

## Methods

![(Supplemental) Figure 1](../clean_data/plots/dnXt_flowchar_v3-2.png)

::: {custom-style="CustomCaption"}
Supplemental Figure 1. Workflow for De novo Transcriptome construction and analysis.
:::

## Generation of PacBio long read RNA sequencing data and Illumina short read RNA sequencing

       Human induced pluripotent stem cells (iPSCs) were differentiated into mature RPE following the culturing protocol in Blenkinsop et al<span class="redoc" id="redoc-citation-17"><span class="redoc" id="redoc-citation-47">[-@blenkinsop_human_2015]</span></span>. RNA was isolated from mature RPE 40 days post differentiation and used for \<illumina library prep\> and \<pacbio library prep\>

## Code availability and software versions.

       To improve reproducibility, we wrote all code used to generate both the data and figures for this paper as Snakemake pipelines. <span class="redoc" id="redoc-citation-18">[-@koster_snakemakescalable_2012]</span> All code (and versions) used for this project is publicly available in the following github repositories: <https://github.com/vinay-swamy/ocular_transcriptomes_pipeline> (main pipeline),

<https://github.com/vinay-swamy/ocular_transcriptomes_longread_analysis> (long read analysis pipeline), <https://github.com/vinay-swamy/ocular_transcriptomes_paper> (figures and tables for this paper), <https://github.com/vinay-swamy/ocular_transcriptomes_shiny> ([ $
$move to David's gh ¶ ]{.comment-start id="2" author="swamyvs" date="2020-06-23T13:48:29Z"}  webapp) . []{.comment-end id="2"}

WE ALSO NEED TO GIVE A GIT TAG FOR EACH TO INDICATE THE "RELEASE" VERSION WITH THIS PAPER (ALSO DEPOSIT IN ZENODO).

## Computational Analyses

All computational analyses performed in this project are run using multiple Snakemake workflows. Each Snakefile contains the exact parameters for all tools and scripts used in each analysis. All Snakefiles are included as supplementary data.(supplementary data files 1-4)

## Analysis of Long Read Data

      PacBio hifi reads were processed into full length, non-chimeric (FLNC) reads using the Pacbio SMRT link software <span class="redoc" id="redoc-citation-19">[ Have to email nisc about this ]{.comment-start id="3" author="swamyvs" date="2020-06-23T13:48:29Z"} [VERSION][]{.comment-end id="3"}. The existing ENCODE longread RNA-seq pipeline (https://github.com/ENCODE-DCC/long-read-rna-pipeline) was rewritten as a Snakemake workflow as follows. Transcripts were aligned to the human genome using minimap2[-@li_minimap2_2018]</span>, using an alignment index built on the gencode v28 primary human genome. Sequencing errors in aligned long reads were corrected using TranscriptClean <span class="redoc" id="redoc-citation-20">[-@wyman_transcriptclean_2019]</span>, using default paramters. Splice junctions for TranscriptClean were obtained using the TranscrtiptClean accessory script "get_SJs_from_gtf.py" using the gencode v28 comprehensive annotation as the input. A list of common variants to avoid correcting were obtained from the ENCODE portal (https://www.encodeproject.org/files/ENCFF911UGW/). The long read transcriptome was generated with TALON <span class="redoc" id="redoc-citation-21">[-@wyman_technology-agnostic_2020]</span>. A TALON database was generated using talon_initialize_database command, with all default parameters, except for the "--5P" and "--3p" parameters. These parameters represent the maximum distance between close 5' start and 3' ends of similar transcript to merge, and were both set to 100 to match parameters used in later tools. Annotation in gtf format was generated using the talon_create_GTF command, and transcript abundance values were generated using the talon_abundance command.

## Analysis of short read RPE data

      We aligned each sample to the Gencode release 28 hg38 assembly using the genomic aligner STAR and sorted the resulting BAM files using samtools sort. [-@frankish_gencode_2019],<span class="redoc" id="redoc-citation-23">[-@dobin_star_2013]</span>,<span class="redoc" id="redoc-citation-24">[-@li_sequence_2009]</span>. For each sorted BAM file, we constructed a per sample base transcriptome using StringTie with the Gencode V28 comprehensive annotation as a guiding annotation [-@frankish_gencode_2019],[-@pertea_stringtie_2015]. All sample transcriptomes were merged with the long read transcriptome using gffcompare<span class="redoc" id="redoc-citation-27"><span class="redoc" id="redoc-citation-29">[-@pertea_gff_2020]</span></span> with default parameters. We note that the default values for the distance to merge similar 5' starts and 3 ends of transcripts in gffcompare is the same to what we chose for TALON. We define the metric construction accuracy, used to evaluate short read transcriptome construction as the following: $Construction\ Accuracy = \frac {short\ read \ transcriptome\ \cap \ long\ read\ transcriptome} {short\ read \ transcriptome}$

## Construction of tissue specific transcriptomes.

       We used studies with healthy, unperturbed RNA-seq samples from <span class="redoc" id="redoc-citation-28"><span class="redoc" id="redoc-inlinecode-2">52</span> distinct subtissue regions of the body, downloaded and performed quality control the pertinent sequencing data from the sequence read archive (SRA) using methods from our previous work[-@swamy_eye_2019]</span>. We constructed a transcriptome for each sample, and merged samples together to create <span class="redoc" id="redoc-inlinecode-3">52</span> subtissue specific transcriptomes. For each tissue-specific transcriptome, we removed transcripts that had an average expression less than 1 Transcripts Per Million (TPM) across all samples of the same tissue type. All tissue specific transcriptomes were merged to form a single unified GTF annotation file to ensure transcript identifiers were the same across tissues. We merged all ocular tissue transcriptomes to generate a separate pan-eye transcriptome.

## Tissue specific transcriptome quantification

       For each resulting tissue specific transcriptome, we extracted transcript sequences using the tool gffread [-@pertea_gff_2020]  and used these sequences to build a tissue-specific quantification index using the index mode of the alignment free quantification tool Salmon. <span class="redoc" id="redoc-citation-30">[-@patro_salmon_2017]</span> For each sample, we quantified transcript expression using the quant mode of salmon, using a sample's respective tissue specific quantification index. We similarly quantified all ocular samples using the pan-eye transcriptome and the Gencode v28 reference transcriptome.

## Annotation of novel exons

       Analysis of novel transcripts was done using a custom Rscript [ This step is part of a larger script might consider separating into its own script. ]{.comment-start id="4" author="swamyvs" date="2020-06-23T13:48:29Z"}  "annotate\_and\_make\_tissue\_gtfs.R"  []{.comment-end id="4"}. First, a comprehensive set of distinct, annotated exons was generated by mergeing exon annotation from gencode, ensembl, UCSC, and refseq. We then defined a novel exon as any exon within our trancriptome that does not excatly match the chromosome, start, end and strand of an annotated exon. Novels exons were classified as followed. Exons were split into 3 catagories: First, last, and middle exons. We extracted all annotated exon start and stop sites from our set of previously annotated exons. Novel middle exons that had an annotated start but an unannotated end were catagorized as a novel alternative 3' end exons and similarly novel middle exons with an unannoated start but annotated end were catagorized as a novel 5' start exons. Novel middle exons whose start and end both matched annotated exon start and ends are considered retained introns. Novel middle exons whose start and end both did not match annotated starts and ends are considered fully novel exons. We then classified novel first and last exons. Novel first exons are first exons whose start is not in the set of annotared exon starts, and novel last exons are terminal exons whose end is not in the set of annotated exon ends. 

## Validation of DNTX with phylop, CAGE data, and polyA signals
       PhyloP <span class="redoc" id="redoc-citation-31"><span class="redoc" id="redoc-citation-57">[-@pollard_detection_2010]</span></span> scores for the phylop 20 way multi species alignment was downloaded from UCSC's FTP server on October 16th, 2019 and converted from bigWig format to bed format using the wig2bed tool in BEDOPs <span class="redoc" id="redoc-citation-32">[-@neph_bedops_2012]</span> . The average score per exon in both the gencode and DNTX annotation was calulated by intersecting exon locations with phylop scores and then averaging the per base score for each exon, using the intersect and groupby tools from the bedtools suite, respectively. Significant difference in mean phylop score was tested with a wilcox rank sum test. 
       CAGE peaks <span class="redoc" id="redoc-citation-33"><span class="redoc" id="redoc-citation-58">[-@noguchi_fantom5_2017]</span></span> were download from the FANTOM FTP server (https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz) on June 15th 2020. Transcriptional start sites (TSS) were extracted from gencode and DNTX annotations; TSS is defined as the start of the first exon of a transcript. Distance to CAGE peaks was calculated using the closest tool in the bedtools suite. Significant difference in mean distance to CAGE peak between DNTX and gencode annotation was tested with a wilcox rank sum test. 
       Polyadenylation signal annotations were downloaded from the polyA site atlas <span class="redoc" id="redoc-citation-34"><span class="redoc" id="redoc-citation-60">[-@herrmann_polyasite_2020]</span></span> (https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz) on June 15th 2020. Transcriptional end sites(TES) were extracted from gencode and DNTX annotations; TES is defined as the end of the terminal exon of a transcript. Distance to polyA signal was calculated using the closest tool in the bedtools <span class="redoc" id="redoc-citation-35">[-@quinlan_bedtools_2010]</span> suite. Significant difference in mean distance to polyA signal was tested with a wilcox rank sum test. 

## Identification of protein coding novel transcripts.

       We identified protein coding transcripts in our unified transcriptome using the TransDecoder suite [-@haas_novo_2013]. We extracted transcript sequences using the util script "gtf_genome_to_cdna_fasta.pl" and used TransDecoder to find a single best open reading frame from each transcript. We then used the "agat_sp_add_start_stop.pl" scripts from the AGAT tool (https://github.com/NBISweden/AGAT/) to identify start and stop codons for each open reading frame. Transcripts with no detectable ORF or missing a start or stop codon were labelled as noncoding. Additionally, novel isoforms whose predicted amino acid sequence's length was either 200 amino acids shorter than the shortest annotated isoform for that gene were marked as noncoding. Similarly, novel isoforms whose predicted amino acid sequence's length was 200 aa greater than the longest annotated isoform were also marked as noncoding.

## Prediction of novel loci function

       Protein coding novel loci were compared to the uniprot protein sequence data base using blastp<span class="redoc" id="redoc-citation-37">[-@altschul_basic_1990]</span>. Blastp results were integrated with hmmer, which identified protein families and domains associated with each novel loci.

## Analysis of novel isoforms in eye tissues.

An Upset <span class="redoc" id="redoc-citation-38">[-@lex_upset_2014]</span> plot was generated using the ComplexUpset package(https://github.com/krassowski/complex-upset). Fraction Isoform Usage(FIU) for each calcualted each transcript t associated with a parent gene g using the following formula: $FIU_t =  \frac{TPM_t}{TPM_g}$ . Raincloud plots of FIU were generated using the R_Rainclouds package. <span class="redoc" id="redoc-citation-39">[-@allen_raincloud_2019]</span>

## Prediction of Variant impact using *de novo* transcriptomes.

Clinvar variants were downloaded from the clinvar database in VCF format and ClinVar metadata in tab delmited format were downloded on April 19th 2020. The VCF of variants was used as the input variants for the Variant Effect Predictor from Ensembl <span class="redoc" id="redoc-citation-40"><span class="redoc" id="redoc-citation-66">[-@mclaren_ensembl_2016]</span></span> . Tissue specific variant priority was determined using gtfs corresponding to each tissue specific transcriptome as input annotation. Variant predictions were also generated using the gencode v28 comprhensive transcript annotation. Variants marked as "Uncertain Significance" in the tab delimited clinvar file were used in down stream analysis. The following query was used to search for variants associated with ocular disease: 'macula|retin|leber|cone|cornea|bardet|ocular|optic|joubert'

## Analysis of fetal retina RNA-seq data.
       RNA-seq samples from Mellough et al. were downloaded from the SRA using methods from a previous study  [-@swamy_eye_2019] . Samples were quantified using salmon with a quantification index generated using our fetal retina *de novo* transcriptome. We removed samples deemed as outliers by first performing principal component analysis of transcript level expression data, calculating the center of all data using the first two principal compenents, and removing the 5 samples furthest away from the center. The remaining samples were normalized using calcNormFactors from the edgeR <span class="redoc" id="redoc-citation-42">[-@robinson_edger_2010]</span> R package and converted to weights using the voom function from the limma  R package. <span class="redoc" id="redoc-citation-43">[-@ritchie_limma_2015]</span> Differential expression was modeled using the lmFit function using developmental time point as the model design and tested for signficant change in expression using the Ebayes function from limma. Gene Set enrichment was testing using the clusterprofileR package. <span class="redoc" id="redoc-citation-44">[-@yu_clusterprofiler_2012]</span> Heatmaps were generated using the ComplexHeatmap package <span class="redoc" id="redoc-citation-45">[-@gu_complex_2016]</span> .

## Computing Resources

       All computation was performed on the National Institutes of Health high performance compute system Biowulf (hpc.nih.gov).

## Figures and Tables

       All statistical analyses, figures and tables in this paper were generated using the R programming language. <span class="redoc" id="redoc-citation-46">[-@r_core_team_r_2019]</span> A full list of packages and versions can be found in supplementary file session\_info.txt

# Results

## Long Read Pacbio RNA sequencing guides *de novo* transcriptome construction

In order to evaluate the accuracy of short read transcriptome construction, we first generated PacBio long read RNA-seq data and Illumina short read RNA-seq data from a stem cell derived RPE cell line. These cell lines were cultured using a highly optimized protocol, and thus should have minimal biological variation [-@blenkinsop_human_2015]. We used this sequencing data to construct a long read transcriptome and a short read transcriptome.(methods) In our long read transcriptome we find <span class="redoc" id="redoc-inlinecode-4">1163239</span> distinct transcripts, and in our short read transcriptome <span class="redoc" id="redoc-inlinecode-5">366888</span> distinct transcripts

<div class="redoc" id="redoc-codechunk-3">


</div>
<div class="redoc" id="redoc-codechunk-4">
![](txome_paper_v2_files/figure-docx/longread_results-1.png)<!-- -->

</div>

::: {custom-style="CustomCaption"}
Figure 1. Transcript lengths are substantially longer in the long-read based approach. A,B) Intersection of novel and annotated transcript lengths between Pacbio (long read) and Stringtie (short read) transcriptomes. The total number of constructed transcripts is given in the text to the right of the violin plot. C) Short read construction accuracy stratified by transcript length at different TPM based transcript exclusion thresholds.
:::

    In our initial comparison between short and and long read transcriptomes, we see a low transcriptome construction accuracy (see Methods) of <span class="redoc" id="redoc-inlinecode-6">0.208</span>. When examining the transcript lengths of each build we see that the two methods show very different transcript length distributions for both novel and previously annotated transcripts, with the short read build comprised mostly of smaller transcripts (Fig 1A). As the PacBio data was generated using two libraries for 2000 bp and >3000 bp (need to double check this), we expected to see enrichment for longer transcripts in the pacbio data set. (Supplemental Fig 2) To assess accuracy relative to transcript length, we group transcripts by length in 1000 bp intervals, and compare accuracy between each group. We found the accuracy significantly improves for transcripts longer than 2000 bp. The construction accuracy is <span class="redoc" id="redoc-inlinecode-7">0.426</span> and <span class="redoc" id="redoc-inlinecode-8">0.137</span> for transcripts above and below 2000 bp, respectively.(Fig 1B)
    We experimented with various methods to remove spurious transcripts to improve construction accuracy. We first removed transcripts that were not expressed at 1 TPM in at least one sample as outlined in the stringtie's recommended protocol. <span class="redoc" id="redoc-citation-48">[-@pertea_transcript-level_2016]</span> This improved construction accuracy to <span class="redoc" id="redoc-inlinecode-9">0.475</span> for transcripts longer than 2000bp and <span class="redoc" id="redoc-inlinecode-10">0.212</span> for transcripts shorter than 2000bp. As this accuracy was still fairly low, we tried a different filtering scheme. We found the simplest approach with high performance was to retain transcripts that had an average TPM above a specific threshold(Fig 1C). In our downstream pipeline we keep transcripts that have at least an average of 1 TPM across all samples of the same subtissue type as this threshold achieved a build accuracy of <span class="redoc" id="redoc-inlinecode-11">0.772</span> for transcripts longer than 2000Bp and retained <span class="redoc" id="redoc-inlinecode-12">48470</span> transcripts within this short read RPE dataset. 

## A rigorous analysis pipeline finds thousands of novel gene isoforms

<div class="redoc" id="redoc-codechunk-5">


</div>
<div class="redoc" id="redoc-codechunk-6">

Tissue   Source    Samples   Studies   Transcriptome Count
-------  -------  --------  --------  --------------------
Retina   Adult         105         8                 49714
RPE      Fetal          49         7                 49967
Cornea   Adult          43         6                 51469
Retina   Fetal          89         6                 66255
RPE      Adult          48         4                 32012
Cornea   Fetal           6         2                 59408

</div>

::: {custom-style="CustomCaption"}
Table 1. Ocular sample dataset overview and transcriptome count. Transcriptome count is defined as the number of unique transcripts expressed in a given tissue type
:::

      We built transcriptomes from <span class="redoc" id="redoc-citation-49"><span class="redoc" id="redoc-inlinecode-13">368</span> published, publicly available ocular tissue RNA-seq samples using an efficient snakemake pipeline. We include both adult and fetal tissue from cornea, retina, and RPE tissues mined from <span class="redoc" id="redoc-inlinecode-14">29</span> different studies (Table 1). Our fetal tissues consist of both human fetal tissues and human induced pluripotent stem cell (iPSC) derived tissue, as stem cell derived tissue has been showed to closely resemble fetal tissue [-@klimanskaya_derivation_2004]</span> . To more accurately determine tissue specificity of novel ocular transcripts, we supplemented our publicly collated normal (non-disease, perturbation) ocular data set with <span class="redoc" id="redoc-citation-50"><span class="redoc" id="redoc-inlinecode-15">877</span> samples across <span class="redoc" id="redoc-inlinecode-16">46</span> body locations from the GTEx project and constructed transcriptomes for each of these body locations [-@gtex_consortium_genetic_2017]</span> . We refer to each distinct body location as a subtissue here after.
      After initial construction of transcriptomes, we found <span class="redoc" id="redoc-citation-51"><span class="redoc" id="redoc-inlinecode-17">183442</span>  previously annotated transcripts  and <span class="redoc" id="redoc-inlinecode-18">6241675</span> novel transcripts detected in at least one of our <span class="redoc" id="redoc-inlinecode-19">1245</span> samples. We define novel as any region of the human genome that has not been previously annotated within the Gencode, Ensembl, UCSC, and Refseq annotation databases. [-@frankish_gencode_2019]</span> , <span class="redoc" id="redoc-citation-52">[-@zerbino_ensembl_2018]</span> , <span class="redoc" id="redoc-citation-53">[-@oleary_reference_2016]</span> After using the filtering methods desribed above, we merged all tissue specific transcriptomes into a single  final transcriptome which contains  <span class="redoc" id="redoc-inlinecode-20">252983</span> distinct transcripts with <span class="redoc" id="redoc-inlinecode-21">87592</span>  previously annotated and <span class="redoc" id="redoc-inlinecode-22">165391</span> novel transcripts, and includes <span class="redoc" id="redoc-inlinecode-23">114.94673</span> megabases of previously unannotated genomic sequence. (Table 1) We refer to the final transcriptome as the DNTX annotation hereafter.
      We split novel transcripts into two categories: novel isoforms, which are novel variations of known genes, and novel loci, which are previously unreported, entirely novel regions of transcribed sequence.(Fig 2B) Novel isoforms are further classified by the novelty of its encoded protein: an isoform with novel open reading frame, a novel isoform with a known ORF, and isoforms with no ORF as noncoding isoforms.(Fig 2A) The number of distinct ORFs is significantly less than the number of transcripts, with <span class="redoc" id="redoc-inlinecode-24">41700</span> previously annotated ORFs and <span class="redoc" id="redoc-inlinecode-25">43130</span> novel ORFs across all tissues. Across all tissues there is an average of <span class="redoc" id="redoc-inlinecode-26">10392.85</span> novel isoforms and <span class="redoc" id="redoc-inlinecode-27">3579.19</span> novel ORFs.

<div class="redoc" id="redoc-codechunk-7">


</div>
<div class="redoc" id="redoc-codechunk-8">
![](txome_paper_v2_files/figure-docx/overall_stats_figure-1.png)<!-- -->

</div>

::: {custom-style="CustomCaption"}
Figure 2. Overview of Novel Isoforms. A. Number of novel gene isoforms, grouped by transcript type. Brain and body represent an average of 13 and 34 distinct subtissues, respectively B. Novel protein coding and noncoding loci. Novel exon composition of novel isoforms, by isoform type labels indicate number of transcripts. C. Classification of novel exon types, stratified by novel isoform type
:::

       Novel isoforms can occur due to an omission of a previously annotated exon, commonly referred as exon skipping or the addition of a unannotated exon which we refer to as a novel exon. We further classify novel exons by the biological process that may be driving their formations: alternative promoter usage driving the addition of novel first exons (NFE) <span class="redoc" id="redoc-citation-54">[-@landry_complex_2003]</span>, alternative polyadenylation driving the addition of novel terminal exons (NTE) <span class="redoc" id="redoc-citation-55">[-@tian_alternative_2017]</span> , and alternative splicing driving the formation of all novel exons that are not the first or last exon <span class="redoc" id="redoc-citation-56">[-@wang_mechanism_2015]</span> . We further classify alternatively spliced exons into their commonly seen patterns, alternative 5' splice site (A5SS), alternative 3' splice site (A3SS), and retained introns (RI). Exons whose entire sequence was unannotated and is not a retained intron are fully novel exons. We note that all three of these mechanisms can lead to exon skipping, so for simplicity we group all novel isoforms resulting from exon skipping together. (considering adding a Supplemental Figure with diagrams of each group to make it easier for ppl) We found that the majority of novel exons with our dataset are novel first and last exons. We see that the majority of A5SS, A3SSs and RI exons lead to novel protein coding isoforms, whereas novel FE/TE more often lead to noncoding isoforms.  (Fig 2C)

## *De novo* transcriptomes match previously published experimental data than existing annotation

      We validated our *de novo* transcriptomes using three independent datasets. First, we evaluated the conservation of our transcriptomes, as conservation has been a historic marker for function. We used PhyloP 20 way species alignment [-@pollard_detection_2010], a measure of conservation between species, to calculate the average conservation score for each exon our DNTX annotation, and compared that to the average conservations score for each exon in the gencode annotation. We found that on average, exons in out DNTX annoation are more conserved than exons in the gencode annotation (pvalue <2.2e-16) (Supplemental Figure 3A). 
      Next, as we saw an enrichment of novel first and last exons within our data set, we decided to compare the transcriptional start sites (TSS) and transcriptional end sites(TES) within our DNTX annotation to two well established annotation databases. We compared DNTX and gencode TSS's to CAGE-seq data from the FANTOM consortium.  [-@noguchi_fantom5_2017]  As CAGE-seq is optimized to detect the 5' end of transcripts, we reason that it can serve as a valid ground truth set to evaluate TSS detection. <span class="redoc" id="redoc-citation-59">[-@takahashi_cage-_2012]</span> We calculated the absolute distance of DNTX TSS's to CAGE peaks, and compared them to the absolute distance of gencode TSS's to CAGE peaks. We found that on average DNTX TSS's are closer to CAGE peaks than gencode TSS's (pvalue <2.2e-16)(Supplemental Figure 3B). Next we evaluated TES's using the polyA atlas, which is comprhensive annotation of polyadenylation signals generated from aggregating 3' seq data from multiple studies. [-@herrmann_polyasite_2020] As 3'-seq data is designed to accurately capture the 3' ends of transcripts it can similarly serve as a ground truth set to evalute the accuracy of TES's. <span class="redoc" id="redoc-citation-61">[-@beck_3-end_2010]</span> We calculated the absolute distance of DNTX TES's to annotated polyA signalsand compared them to the absolute distance of gencode TES's to polyA signals. We found that on average DNTX TES's are closer to annotated polyadenylation signals than gencode TSS's (pvalue <2.2e-16) (Supplemental Figure 3C)

## *De novo* transcriptomes reduce overall transcriptome sizes

      Our transcriptomes removed on average <span class="redoc" id="redoc-inlinecode-28">76.141</span> % of a tissue's base transcriptome. We define base transcriptome for a tissues as any transcript in the gencode annotation with non zero TPM in at least one sample of a given tissue type. This was a large reduction in transcriptome size and ee wanted to ensure we were not unduly throwing away data. We quantified transcript expression of our samples using Salmon, quantifying each sample twice: once using the full gencode V28 human transcript annotation, and once using its associated tissue specific transcriptome We found that despite the <span class="redoc" id="redoc-inlinecode-29">76.141</span> % reduction in number of transcripts between the base gencode and *de novo* transcriptomes (Supplemental Figure 4A), the average salmon mapping rate increases  by  <span class="redoc" id="redoc-inlinecode-30">2.041</span> %  indicating that the vast majority of gene expression data is retained within our transcriptome. (Supplimental Figure 4B)

## Novel Isoforms in Ocular tissues

<div class="redoc" id="redoc-codechunk-9">


</div>
<div class="redoc" id="redoc-codechunk-10">
![](txome_paper_v2_files/figure-docx/novel_isoforms_ocular_tissues-1.png)<!-- -->

</div>

::: {custom-style="CustomCaption"}
Figure 3. Analysis of novel isoform characteristics A). Set intersection of novel isoforms in ocular transcriptomes. B). Boxplots of fraction isoform usage(FIU) overlaid over FIU data points with estimated distribution of data set above each boxplot
:::

      Next, we analyzed the novel isoforms within our ocular transcriptomes. We compared the overlap in constructed novel isoforms across ocular tissues and found that <span class="redoc" id="redoc-inlinecode-31">94.99</span> % of novel isoforms are specific to a singular ocular subtissue (Fig 3A). For each novel isoform we then calculated fraction isoform usage (FIU), or the fraction of total gene expression a transcript contributed to its parent gene. We found that on average novel isoforms contribute to <span class="redoc" id="redoc-inlinecode-32">25.79</span> % of their parent gene's expression.

## Differential Usage of Gene Isoforms Occurs during Retinal Development

Mellough et al showed that alternative splicing plays a role during retinal development <span class="redoc" id="redoc-citation-62">[-@mellough_integrated_2019]</span> As this is a subset of alternative isoform usage, we hypothesized that alternative isoform usage plays a role in retinal development. We use RNA-seq data from  Mellough et al that we did not include in the data used to build our transcriptomes. We used our fetal retina *de novo* transciptome to quantify transcript expression and analyzed differential transcript usage (DTU) across development.

<div class="redoc" id="redoc-codechunk-11">
![](txome_paper_v2_files/figure-docx/fetal_retina_diff_exp-1.png)<!-- -->

</div>

::: {custom-style="CustomCaption"}
Figure 4 Differential Transcript usage during Retinal Development. A) Volcano Plot of tested transcripts B) Dot plot for gene set enrichment analysis E) Heatmap of genes with DTU associated with eye development  D) Transcript models for MYO9A, a gene undergoing DTU F) change in MYO9A FIU across development F) average TPM expression of MYO9A across development 
:::

We analyzed <span class="redoc" id="redoc-citation-63"><span class="redoc" id="redoc-inlinecode-33">24</span> samples across <span class="redoc" id="redoc-inlinecode-34">14</span> developmental days post fertilization and found <span class="redoc" id="redoc-inlinecode-35">1717</span> transcripts acrosss <span class="redoc" id="redoc-inlinecode-36">812</span> genes involved in diffferential transcript usage(DTU).(Fig 4A) We define DTU as an transcript that is differentially expressed (qvalue <.01) and has a FIU difference of .25 in at least one comparison of time points. We found that genes involved in DTU were enriched(qvalue <.05) for genes related to eye and neurological development.(Fig 4B), and that hiearchical clustering of DTU transcripts geneates an early stage and late stage cluster.(Fig 4C) One of these genes, MYO9A, is a perfect example of DTU. MYO9A is associated with the visual perception GO term and plays a role in ocular development and has been associated with ocular disease [-@gorman_cloning_1999]</span> . While expression of MYO9A remains relatively unchanged across development, expression of two of its associated isoforms(Fig 2D) changes across during development: a novel isoforms is highly expressed early during development, but switches to a cannonical isoform later in development.(Fig 2E,F)

## *de novo* transcriptome allow for a more precise variant prioritization.

The prediction of variant's biological impact is a fundamental step in variant prioritization, which is a critical step in the diagnosis of genetic disease and in the interpretation GWAS results. The recent GnomAD<span class="redoc" id="redoc-citation-64">[-@karczewski_mutational_2020]</span> project showed that taking transcript expression in to account can improve variant prioritization. <span class="redoc" id="redoc-citation-65">[-@cummings_transcript_2020]</span> We hypothesize that our *de novo* transcriptomes can provide an additional layer of refinement to expression based variant conseqeunce prediction. We use Ensembl's Variant Effect Predictor combined with our tissue specific transcriptomes to generate a tissue specific prediction of variant impact, and reclassify variants of unknown significance (VUS) in the ClinVar database. [-@mclaren_ensembl_2016],<span class="redoc" id="redoc-citation-67">[-@landrum_clinvar_2018]</span>. We grouped different predicted consequences into a High, moderate and low impact groups and provide a full mapping of biological consequence to impact in our supplementary data.

<div class="redoc" id="redoc-codechunk-12">


</div>
<div class="redoc" id="redoc-codechunk-13">
![](txome_paper_v2_files/figure-docx/vep-1.png)<!-- -->

</div>

::: {custom-style="CustomCaption"}
Figure 5. Tissue specific variant prioritization of two variants of unknown significance (VUS) previously associated with ocular disease. A) Variant Impact dot plot. Size corresponds to average TPM of the transcript with the highest TPM  and color to impact. A full mapping of variant descriptors to impact can be found in supplemental data. b) Transcript Models of CACNA2D4 thick lines represent CDS; Asterisk represents location of variant 631715. Arrows indicate direction of transcription
:::

We compared our tissue specific variant effect predictions to variant effect predictions generated using the gencode transcript annotation. We found that globally out of <span class="redoc" id="redoc-inlinecode-37">226824</span> total variants examined, the priority of <span class="redoc" id="redoc-inlinecode-38">118297</span> variants decreased in at least one tissue, and <span class="redoc" id="redoc-inlinecode-39">1324</span> variants increased in at least one tissue. We see a larger amount of variants decrease in predicted effects as the transcript which previously lead to an effect is not expressed in a particular tissue. To demonstrate the utility of tissue specific variant prioitizaion, we focused specifically on variants previously associated with ocular disease. We examined ClinVar variant 631672, which was previously associated with retinal cone dystrophy 4, and has a high predicted impact only in ocular tissues, and not in other tissues. (Fig 5A) Variant 631672 is in the gene CACNA2D4 and in all transcripts expressed in ocular tissues this variant leads to a premature stop codon. While CACNA2D4 has other isoforms expressed in other tissues, these isoforms are not affected by variant 631672 and thus there is a higher probability that this variant may be playing a role in disease.(Fig 5B) Conversely, in clinvar variant 631715 which is associated with posterior column ataxia-retinitis pigmentosa syndrome and is within gene FLVCR2, the variant is low impact on transcripts expressed in eye tissues, but higher impact in other body tissues. This may suggest that this variant is less likely to be associated with the disease.

## A companion visualization tool enables easy use of *de novo* transcriptomes

<div class="redoc" id="redoc-codechunk-14">
![](txome_paper_v2_files/figure-docx/app_viz-1.png)<!-- -->

</div>

::: {custom-style="CustomCaption"}
Figure 6. Screenshots from dynamic *de novo* transcriptome visualization tool. A). FIU bar plot for selected gene and tissue. B). Exon level diagram of transcript body Thicklines represent coding region of transcript. novel exons colored in red. Tooltip contains genomic location and phylop score C) Bargraph of fraction of samples within dataset each transcript was consructed in by tissue.
:::

       To make our results easily accessible, we designed a web-app for visualizing and accessing our *de novo* transcriptomes. Users start by selecting a gene or searching for a gene by genomic location, and can choose up to 5 tissues to visualize transcript expression in. For each tissue we show the FIU for each transcript associated with a gene. (Fig 6A) We show the exon-intron structure of each transcript and mousing over exons show genomic location overlapping SNPs, and phylogentic conservation score. (Fig 6B) We additionally show a barplot of the fraction of samples in each tissue each transcript was constructed in.(Fig 6C) Users can  download the *de novo* transcriptomes for selected tissues in GTF and fasta format. While visualization of direct transcript expresion is not a part of this app, it can be viewed in the eyeIntegration app [-@swamy_eye_2019] by selected 'DNTX' as the transcript annotation.

## A easy-to-use pipeline to enable future research

[ . (havent worked it out entirely, but snakemake can automatically use a singulaity container, so if I provide a singularity image anyone can use it) ]{.comment-start id="5" author="swamyvs" date="2020-06-23T13:48:29Z"}  Finally, we package all tools used for our transcriptome pipeline within a portable docker container with a stand-alone run script. This pipeline allows other researchers to run their own samples, and generate figures and annotations similar to what is shown here []{.comment-end id="5"}

# Discussion

      [ DISCUSSION IS TOO LONG CUT DOWN ]{.comment-start id="6" author="swamyvs" date="2020-06-23T13:48:29Z"} Motivated  []{.comment-end id="6"} by the lack of a comprehensive pan eye transcriptome we created the first comprehensive set of ocular transcriptomes. By using long read RNA-sequencing data to calibrate our short-read construction pipeline, we increase our confidence that we are creating real, biologically relevant transcriptomes. We found that concordance between long and short read based transcriptome is directly related to transcript length and transcript expression across many tissues. We saw a clear inability within this PacBio data set to accurately detect transcripts shorter than 2000Bp for both previously annotated and novel transcripts. As many of the transcripts constructed using short reads are below this threshold, long read sequencing data enriched for smaller transcript sizes would provide greater insight in future studies. 

      We used a large dataset compiled from published RNA-seq data to build our ocular transcriptomes, an approach which has several key advantages. First, our large sample size allows us to combat the noisy nature of RNA-seq data. Second, as our samples come from mutliple studies and types of sample preparation because , we can be more confident that our transcriptomes accurately reflect the biology of its originating subtissue and are not a technical artifact due to preparation of the samples. Further more, we are confident that these results are accurate because our *de novo* transcritomes match existing large scale data sets and are more conserved than existing annotation. (Supplemental Figure 2) 

      In each tissue we examined, we found hundereds of novel gene isoforms, many of which were novel due to novel exons. Within ocular tissues, these novel isoforms are most commonly specific to single subtissue. This makes sense as over half of the exons in our *de novo* transcirptomes are first and last exons, which have been prevously shown to significanlty contribute to the tissue specificity of gene isoforms <span class="redoc" id="redoc-citation-69">[-@reyes_alternative_2018]</span> We also find that on average novel isoforms represent about <span class="redoc" id="redoc-citation-70"><span class="redoc" id="redoc-inlinecode-40">25.79</span> % of their parent gene's expression. However, it is difficult to say if these are non-functional. Its is entirely possible that these isoforms are from low population cell types, as transcript annotation was previously shown to be incopmlete in rare cell types. [-@zhang_incomplete_2020]</span> This especially makes sense in the retina which contain multiple distinct cell types, several of which contribute to 5% or less to total cell population <span class="redoc" id="redoc-citation-71">[-@yan_cell_2020]</span> However, as we imposed a strict expression filter as part of out transcriptome pipeline, we may be removing transcripts specific to rare cell types.

      {== We showed that our ocular transcriptomes are useful in understanding the role of gene isoforms in ocular biology and disease. We saw that multiple novel isoforms of genes associated with retinal development that exhibit significant changes in trascript expression and FIU during retinal development. We used our tissue specific transcriptomes to re prioritize variants of unknown significance, and found the tissue specific variant prioritization can increase the resolution of tissue specific variant prioritization. ==} [ I feel like this needs more but I cant really figure our what to say ]{.comment-start id="7" author="swamyvs" date="2020-06-23T13:48:29Z"}[]{.comment-end id="7"} 

      This work is most useful as a starting point for other researchers; We want to make our transcriptomes easily accessible to other researchers, so we designed a webapp to visualize our transcriptomes and access tissue-specific annotation files. We wanted to provide as much information to the user so that they can draw their own conclusions about the significance of potential novel exon, and so we provide the gene model with novel exon, coding and noncoding regions marked, along with the FIU for each transcript constructed within a gene. We also provide the fraction of samples within a given subtissue type a sample was detected in, to provide a further level of evidence to the validity of constructed transcripts.


# Supplemental Figures

<div class="redoc" id="redoc-codechunk-15">
![](txome_paper_v2_files/figure-docx/sup_fig1 -1.png)<!-- -->

</div>

::: {custom-style="CustomCaption"}
Suppplemental Figure 2. Distribution of PacBio long ead lengths for two library sizes
:::

<div class="redoc" id="redoc-codechunk-16">
![](txome_paper_v2_files/figure-docx/sup_fig2-1.png)<!-- -->

</div>

::: {custom-style="CustomCaption"}
Suppplemental Figure 3. Comparison of DNTX annotation to Gencode Annotation. A) Per exon Phlop score for gencode and DNTX transcripts. B) Average distance of DNTX Transcriptional Start Sites (TSS) and Gencode TSS to CAGE-seq peaks from the FANTOM consortium. C) Average distance of DNTX Transcriptional End Sites (TES) and Gencode TES to polyadenylation signals in the PolyA site atlas.
:::

<div class="redoc" id="redoc-codechunk-17">


</div>
<div class="redoc" id="redoc-codechunk-18">
![](txome_paper_v2_files/figure-docx/map_rate_diff-1.png)<!-- -->

</div>

::: {custom-style="CustomCaption"}
Suppplemental Figure 4. Comparison of Salmon mapping rate change  vs transcriptome size decrease. 
:::

<div class="redoc" id="redoc-codechunk-19">



</div>



# References
