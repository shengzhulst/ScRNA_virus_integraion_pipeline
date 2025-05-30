import glob
import os

# you need to have the "ScRNA_HPV_integration_detection" conda environment installed use the YML file

##########
###set up basic environment parameters
###############
# human virus fasta and gtf files
human_fasta=''
human_gtf=''
human_plus_virus_fasta=''
virus_fasta=''
# this is a file listed all your virus length, an example is here:
# virus  start_pos length 
# HPV16    0   7905
# HPV18    0   7857
# HPV31    0   7912
virus_bed_file=''
## the directory where you put the folder nf-viral-starsolo
dir_for_scripts='./'
#where you put the generated reference files: specifically human and human plus virus
human_reference=''
human_virus_reference=''
#where you put the fastq files, and output files
fq1=''
fq2=''
output_dir=''
# set up your sample ID, varies depending on your data
samples_interest=''
# set up the parameters for the pipeline(notice: only basic parameters are set here, you can change them in the rules)
# number of threads to use
threads=
#length of the reads, usually 100 or 150
read_length=
#length of the UMI, usually 10 or 12
soloUMI_length=
# set the start position for the cell barcode
soloC_start=1
# set the length of the cell barcode, here is an example 16
soloCB_length=16
# set the start position for the UMI
soloUMI_start=17
#whitelist file for cellranger
whitelist='./cellranger-7.2.0/lib/python/cellranger/barcodes/737K-august-2016.txt'


#######################
#### rule part
########################


rule all:
    input:
        # build_reference_HPV_human
        expand(f"{human_reference}"),
        expand(f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.hostAligned.sortedByCoord.out.bam",samples_interest=samples_interest),
        expand(f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.plusAligned.sortedByCoord.out.bam",samples_interest=samples_interest),
        expand(f"{output_dir}"+"star_all/{samples_interest}/human_virus_chimeric_junction.txt",samples_interest=samples_interest),
        expand(f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.validate_inserts.Aligned.sortedByCoord.out.bam",samples_interest=samples_interest),
        expand(f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.vif.refined.tsv",samples_interest=samples_interest),
        expand(f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.VirusDetect.init.genome_plot.png",samples_interest=samples_interest),

        

rule build_reference_HPV_human:
    input:
        human_virus_fasta=f'{human_plus_virus_fasta}',
        human_fasta=f'{human_fasta}'
    output:
        output_reference_human=directory(f'{human_reference}'),
        output_reference_human_virus=directory(f'{human_virus_reference}')
    conda:
        "ScRNA_HPV_integration_detection"
    params:
        gtf_file=f'{human_gtf}',
        sjdboverhang=f'{read_length}-1',
        threads=f'{threads}'
    shell:
        """
        STAR \\
        --runMode genomeGenerate \\
        --genomeDir {output.output_reference_human} \\
        --genomeFastaFiles {input.human_fasta} \\
        --sjdbGTFfile {params.gtf_file} \\
        --runThreadN {params.threads} \\
        --genomeSAindexNbases 14\\
        --limitSjdbInsertNsj 10000000 \\
        --limitGenomeGenerateRAM 77209411328 \\
        --sjdbOverhang {params.sjdboverhang} && \\
        STAR \\
        --runMode genomeGenerate \\
        --genomeDir {output.output_reference_human_virus} \\
        --genomeFastaFiles {input.human_virus_fasta} \\
        --sjdbGTFfile {params.gtf_file} \\
        --runThreadN {params.threads} \\
        --genomeSAindexNbases 14\\
        --limitSjdbInsertNsj 10000000 \\
        --limitGenomeGenerateRAM 77209411328 \\
        --sjdbOverhang {params.sjdboverhang} 
        """

rule run_star_solo_human:
    input:
        fq1=f"{fq1}",
        fq2=f"{fq2}"
    output:
        output_file=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.hostAligned.sortedByCoord.out.bam"
    conda:
        "ScRNA_HPV_integration_detection"
    params:
        human_reference=f"{human_reference}",
        gtf_file=f"{human_gtf}",
        threads=f"{threads}",
        tmp_dir=f"{output_dir}/"+"temp_files/{samples_interest}_tmp/",
        whitelist=f"{whitelist}",
        output_dir=f"{output_dir}"+"star_all/{samples_interest}/",
        sample="{samples_interest}",
        soloUMIlen=f"{soloUMI_length}"
        soloCBstart=f"{soloC_start}",
        soloCBlen=f"{soloCB_length}",
        soloUMIstart=f"{soloUMI_start}"
    shell: '''STAR \
    --genomeDir {params.human_reference} \
    --readFilesIn {input.fq2} {input.fq1}  \
    --outTmpDir {params.tmp_dir} \
    --runThreadN {params.threads} \
    --outFileNamePrefix {params.output_dir}{params.sample}.host \
    --sjdbGTFfile {params.gtf_file} \
    --outSAMattrRGline ID:{params.sample}.host 'SM:{params.sample}.host' 'PL:illumina'  \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outSAMunmapped Within \
    --twopassMode Basic \
    --alignSJDBoverhangMin 5 \
    --genomeSuffixLengthMax 10000 \
    --limitBAMsortRAM 47271261705 \
    --alignInsertionFlush Right \
    --alignIntronMax 100000 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --limitSjdbInsertNsj 2000000 \
    --outReadsUnmapped Fastx \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist {params.whitelist} \
    --soloCBstart {params.soloCBstart} \
    --soloCBlen {params.soloCBlen} \
    --soloUMIstart {params.soloUMIstart} \
    --soloUMIlen {params.soloUMIlen} \
    --clipAdapterType CellRanger4 \
    --soloFeatures Gene GeneFull SJ Velocyto \
    --soloCellFilter EmptyDrops_CR \
    --soloMultiMappers EM \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \

if [ -f {params.output_dir}{params.sample}.hostUnmapped.out.mate1 ]; then
    mv {params.output_dir}{params.sample}.hostUnmapped.out.mate1 {params.output_dir}{params.sample}.host.unmapped_1.fastq
    gzip -f {params.output_dir}{params.sample}.host.unmapped_1.fastq
fi &&
if [ -f {params.output_dir}{params.sample}.hostUnmapped.out.mate2 ]; then
    mv {params.output_dir}{params.sample}.hostUnmapped.out.mate2 {params.output_dir}{params.sample}.host.unmapped_2.fastq
    gzip -f {params.output_dir}{params.sample}.host.unmapped_2.fastq
fi   
'''

rule run_star_solo_human_HPV:
    input:
        fq1=f"{output_dir}"+"star_all/{samples_interest}/"+"{samples_interest}.host.unmapped_1.fastq.gz",
        fq2=f"{output_dir}"+"star_all/{samples_interest}/"+"{samples_interest}.host.unmapped_2.fastq.gz",
    output:
        output_file=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.plusAligned.sortedByCoord.out.bam"
    conda:
        "ScRNA_HPV_integration_detection"
    params:
        human_virus_reference=f"{human_virus_reference}",
        gtf_file=f"{human_gtf}",
        threads=f"{threads}",
        tmp_dir=f"{output_dir}"+"temp_files/{samples_interest}_virus_tmp/",
        whitelist=f"{whitelist}",
        output_dir=f"{output_dir}"+"star_all/{samples_interest}/",
        sample="{samples_interest}",
        soloUMIlen=f"{soloUMI_length}",
        soloCBstart=f"{soloC_start}",
        soloCBlen=f"{soloCB_length}",
        soloUMIstart=f"{soloUMI_start}"
    shell:"""STAR \
    --genomeDir {params.human_virus_reference} \
    --readFilesIn {input.fq1} {input.fq2}  \
    --outTmpDir {params.tmp_dir} \
    --runThreadN {params.threads} \
    --outFileNamePrefix {params.output_dir}{params.sample}.plus \
    --sjdbGTFfile {params.gtf_file} \
    --outSAMattrRGline ID:{params.sample}.plus 'SM:{params.sample}.plus' 'PL:illumina'  \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outSAMunmapped Within \
    --twopassMode Basic \
    --alignSJDBoverhangMin 5 \
    --genomeSuffixLengthMax 10000 \
    --limitBAMsortRAM 47271261705 \
    --alignInsertionFlush Right \
    --alignIntronMax 100000 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --limitSjdbInsertNsj 2000000 \
    --chimJunctionOverhangMin 12 \
    --chimOutJunctionFormat 0 \
    --chimSegmentMin 10 \
    --chimSegmentReadGapMax 3 \
    --chimScoreJunctionNonGTAG 0 \
    --chimNonchimScoreDropMin 10 \
    --chimMultimapScoreRange 10 \
    --chimMultimapNmax 2 \
    --chimOutType Junctions WithinBAM \
    --outReadsUnmapped Fastx \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist {params.whitelist} \
    --soloCBstart {params.soloCBstart} \
    --soloCBlen {params.soloCBlen} \
    --soloUMIstart {params.soloUMIstart} \
    --soloUMIlen {params.soloUMIlen} \
    --clipAdapterType CellRanger4 \
    --soloFeatures Gene GeneFull SJ Velocyto \
    --soloCellFilter EmptyDrops_CR \
    --soloMultiMappers EM \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \

if [ -f {params.output_dir}{params.sample}.plusUnmapped.out.mate1 ]; then
    mv {params.output_dir}{params.sample}.plusUnmapped.out.mate1 {params.output_dir}{params.sample}.plus.unmapped_1.fastq
    gzip -f {params.output_dir}{params.sample}.plus.unmapped_1.fastq
fi
if [ -f {params.output_dir}{params.sample}.plusUnmapped.out.mate2 ]; then
    mv {params.output_dir}{params.sample}.plusUnmapped.out.mate2 {params.output_dir}{params.sample}.plus.unmapped_2.fastq
    gzip -f {params.output_dir}{params.sample}.plus.unmapped_2.fastq
fi
samtools index -@ 8 {params.output_dir}{params.sample}.plusAligned.sortedByCoord.out.bam
        """

rule run_python_script_step1:
    input:
        junction=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.plusChimeric.out.junction",
        bam_plus=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.plusAligned.sortedByCoord.out.bam",
        bam_host=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.hostAligned.sortedByCoord.out.bam"

    output:
        chimj=f"{output_dir}"+"star_all/{samples_interest}/human_virus_chimeric_junction.txt",

    conda:
        "ScRNA_HPV_integration_detection"
    params:
        fasta_HPV=f'{virus_fasta}',
        output_dir=f"{output_dir}"+"star_all/{samples_interest}/",
        sample="{samples_interest}",
        human_fasta=f'{human_fasta}',
        dir_for_scripts=f'{dir_for_scripts}'
    shell:
        """
        sed '1 s/$/\tbarcode\tUMI/' {input.junction} > {params.output_dir}{params.sample}.plus.addheader.Chimeric.out.junction &&
        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/pre_filter_non_human_virus_chimeric_alignments.py --chimJ {params.output_dir}{params.sample}.plus.addheader.Chimeric.out.junction \
        --viral_db_fasta {params.fasta_HPV} \
        --output {output.chimj} &&
        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/chimJ_to_virus_insertion_candidate_sites.py \
        --chimJ {output.chimj} \
        --viral_db_fasta {params.fasta_HPV} \
        --max_multi_read_alignments 10 \
        --output_prefix {params.output_dir}{params.sample}.vif.init.tmp \
        --remove_duplicates &&
        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/extract_prelim_chimeric_genome_read_alignments.py \
        --star_bam {params.output_dir}{params.sample}.plusAligned.sortedByCoord.out.bam \
        --vif_full_tsv {params.output_dir}{params.sample}.vif.init.tmp.full.tsv \
        --output_bam {params.output_dir}{params.sample}.vif.init.genome_chimeric_evidence.bam &&
        python3 {params.dir_for_scripts}nf-viral-starsolo/greedily_assign_multimapping_reads_among_insertions.py \
        --init_full_tsv {params.output_dir}{params.sample}.vif.init.tmp.full.tsv \
        --include_readnames > {params.output_dir}{params.sample}.vif.init.tmp.full.virus_breakend_grouped.tsv &&
        python3 {params.dir_for_scripts}nf-viral-starsolo/viralintegration/bin/incorporate_read_alignment_stats.py \
        --supp_reads_bam {params.output_dir}{params.sample}.vif.init.genome_chimeric_evidence.bam \
        --vif_full_tsv {params.output_dir}{params.sample}.vif.init.tmp.full.virus_breakend_grouped.tsv \
        --output {params.output_dir}{params.sample}.vif.init.full.w_read_stats.tsv &&
        python3 {params.dir_for_scripts}nf-viral-starsolo/modified_files_from_viralintegration-0.1.1/bin/incorporate_breakpoint_entropy_n_splice_info.py \
        --vif_tsv  {params.output_dir}{params.sample}.vif.init.full.w_read_stats.tsv \
        --ref_genome_fasta {params.human_fasta} \
        --viral_genome_fasta {params.fasta_HPV} \
        --output {params.output_dir}{params.sample}.vif.init.full.w_brkpt_seq_entropy.tsv &&
        python3 {params.dir_for_scripts}nf-viral-starsolo/viralintegration/bin/revise_primary_target_list_via_brkpt_homologies.py \
        --vif_tsv {params.output_dir}{params.sample}.vif.init.full.w_brkpt_seq_entropy.tsv > {params.output_dir}{params.sample}.vif.init.full.tsv
        ########rm {params.output_dir}{params.sample}.vif.init.tmp.full.tsv &&
        python3 {params.dir_for_scripts}nf-viral-starsolo/abridge.py \
        -i {params.output_dir}{params.sample}.vif.init.full.tsv \
        -o {params.output_dir}{params.sample}.vif.init.filtered.tsv \
        -a {params.output_dir}{params.sample}.vif.init.filtered.abridged.tsv && 
        python3 {params.dir_for_scripts}nf-viral-starsolo/viralintegration/bin/extract_chimeric_genomic_targets.py \
        --fasta {params.human_fasta} \
        --patch_db_fasta {params.fasta_HPV} \
        --output_prefix {params.output_dir}{params.sample}.vif.extract \
        --chim_events {params.output_dir}{params.sample}.vif.init.filtered.abridged.tsv \
        --pad_region_length 500
        """




rule run_star_solo_human_back_again:
    input:
        fq1=f"{output_dir}"+"star_all/{samples_interest}/"+"{samples_interest}.host.unmapped_1.fastq.gz",
        fq2=f"{output_dir}"+"star_all/{samples_interest}/"+"{samples_interest}.host.unmapped_2.fastq.gz",
    output:
        output_file=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.validate_inserts.Aligned.sortedByCoord.out.bam"
    conda:
        "ScRNA_HPV_integration_detection"
    params:
        human_virus_reference=f"{human_virus_reference}",
        human_reference=f"{human_reference}",
        gtf_file=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.vif.extract.gtf",
        threads=14,
        #tmp_dir=f"{output_dir}"+"{samples_interest}_back_human_tmp/",
        tmp_dir=f"{output_dir}"+"temp_files/{samples_interest}_back_human_tmp/",
        whitelist=f"{whitelist}",
        output_dir=f"{output_dir}"+"star_all/{samples_interest}/",
        sample="{samples_interest}",
        fasta_used=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.vif.extract.fasta",
        soloUMIlen=f"{soloUMI_length}",
        soloCBstart=f"{soloC_start}",
        soloCBlen=f"{soloCB_length}",
        soloUMIstart=f"{soloUMI_start}"

    shell:
        """
    STAR \
    --genomeDir {params.human_reference} \
    --readFilesIn {input.fq1} {input.fq2}  \
    --outTmpDir {params.tmp_dir} \
    --runThreadN {params.threads} \
    --outFileNamePrefix {params.output_dir}{params.sample}.validate_inserts.  \
    --genomeFastaFiles {params.fasta_used} \
    --outSAMattrRGline ID:{params.sample}.validate_inserts 'SM:{params.sample}.validate_inserts.' 'PL:illumina'  \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outSAMunmapped Within \
    --twopassMode Basic \
    --alignSJDBoverhangMin 5 \
    --genomeSuffixLengthMax 10000 \
    --limitBAMsortRAM 47271261705 \
    --alignInsertionFlush Right \
    --outSAMfilter KeepOnlyAddedReferences \
    --alignMatesGapMax 100000 \
    --alignIntronMax 100000 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --limitSjdbInsertNsj 2000000 \
    --outReadsUnmapped Fastx \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist {params.whitelist} \
    --soloCBstart {params.soloCBstart} \
    --soloCBlen {params.soloCBlen} \
    --soloUMIstart {params.soloUMIstart} \
    --soloUMIlen {params.soloUMIlen} \
    --clipAdapterType CellRanger4 \
    --soloFeatures Gene GeneFull SJ Velocyto \
    --soloCellFilter EmptyDrops_CR \
    --soloMultiMappers EM \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR 

if [ -f {params.output_dir}{params.sample}.validate_inserts.Unmapped.out.mate1 ]; then
    mv {params.output_dir}{params.sample}.validate_inserts.Unmapped.out.mate1 {params.output_dir}{params.sample}.validate_inserts.unmapped_1.fastq
    gzip -f {params.output_dir}{params.sample}.validate_inserts.unmapped_1.fastq
fi
if [ -f {params.output_dir}{params.sample}.validate_inserts.Unmapped.out.mate2 ]; then
    mv {params.output_dir}{params.sample}.validate_inserts.Unmapped.out.mate2 {params.output_dir}{params.sample}.validate_inserts.unmapped_2.fastq
    gzip -f {params.output_dir}{params.sample}.validate_inserts.unmapped_2.fastq
fi
samtools index -@ 8 {params.output_dir}{params.sample}.validate_inserts.Aligned.sortedByCoord.out.bam
        """


rule run_python_script_step2:
    input:
        bam_plus=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.validate_inserts.Aligned.sortedByCoord.out.bam"

    output:
        refine_tsv=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.vif.refined.tsv"
    conda:
        "ScRNA_HPV_integration_detection"
    params:
        fasta_HPV=f'{virus_fasta}',
        output_dir=f"{output_dir}"+"star_all/{samples_interest}/",
        sample="{samples_interest}",
        human_fasta=f'{human_fasta}',
        dir_for_scripts=f'{dir_for_scripts}',
        gtf_file=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.vif.extract.gtf",
        virus_bed=f"{virus_bed_file}",
        dedup_bam=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.dedup.bam"
    shell:
        """
        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/bam_mark_duplicates.py \
        -i {input.bam_plus} \
        -o {params.dedup_bam} \
        -r &&
        samtools index {params.dedup_bam} &&
        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/chimeric_contig_evidence_analyzer.py \
        --patch_db_bam {params.dedup_bam} \
        --patch_db_gtf {params.gtf_file} \
        --output_prefix {params.output_dir}{params.sample}.vif &&
        samtools index {params.output_dir}{params.sample}.vif.evidence.bam &&
        samtools view -b -L {params.virus_bed} -o {params.output_dir}{params.sample}.VirusDetect.igvjs.bam {params.output_dir}{params.sample}.plusAligned.sortedByCoord.out.bam &&
        samtools sort -o {params.output_dir}{params.sample}.VirusDetect.sorted.igvjs.bam -T {params.output_dir}{params.sample}.VirusDetect {params.output_dir}{params.sample}.VirusDetect.igvjs.bam &&
        samtools index {params.output_dir}{params.sample}.VirusDetect.sorted.igvjs.bam &&
        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/restrict_bam_to_proper_aligns.py {params.output_dir}{params.sample}.VirusDetect.sorted.igvjs.bam {params.output_dir}{params.sample}.VirusDetect.sorted.igvjs.bam.clean.bam &&
        mv {params.output_dir}{params.sample}.VirusDetect.sorted.igvjs.bam.clean.bam {params.output_dir}{params.sample}.VirusDetect.sorted.igvjs.bam &&
        mv {params.output_dir}{params.sample}.VirusDetect.sorted.igvjs.bam.clean.bam.bai {params.output_dir}{params.sample}.VirusDetect.sorted.igvjs.bam.bai &&
        Rscript {params.dir_for_scripts}/nf-viral-starsolo/modified_files_from_viralintegration-0.1.1/bin/refine_VIF_output.R --prelim_counts {params.output_dir}{params.sample}.vif.init.filtered.abridged.tsv \
        --vif_counts {params.output_dir}{params.sample}.vif.evidence_counts.tsv \
        --output {params.output_dir}{params.sample}.vif.prelim.refined.tsv &&
        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/examine_flanking_uniq_kmer_composition.py \
        --vif_tsv {params.output_dir}{params.sample}.vif.prelim.refined.tsv \
        --min_frac_uniq 0.0 \
        --output {params.output_dir}{params.sample}.vif.refined.tsv &&

        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/distill_to_primary_target_list_via_brkpt_homologies.py \
        --vif_tsv {params.output_dir}{params.sample}.vif.refined.tsv > {params.output_dir}{params.sample}.vif.refined.distilled.tsv
        """


rule run_visulization_script_part:
    input:
        input_tsv=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.vif.init.filtered.abridged.tsv"
    
    output:
        output_png=f"{output_dir}"+"star_all/{samples_interest}/{samples_interest}.VirusDetect.init.genome_plot.png"
    
    conda:
        "ScRNA_HPV_integration_detection"
    params:
        output_dir=f"{output_dir}"+"star_all/{samples_interest}/",
        sample="{samples_interest}",
        virus_fai=f"{virus_fasta}.fai",
        fasta_HPV=f"{virus_fasta}",
        gtf_file=f"{human_gtf}",
    shell:
        """    
        Rscript {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/make_VIF_genome_abundance_plot.R \
        --vif_report {input.input_tsv} \
        --output_png {params.output_dir}{params.sample}.VirusDetect.init.genome_plot.png \
        --title 'Preliminary_Genome_Wide_Abundance' &&
        Rscript {params.dir_for_scripts}/nf-viral-starsolo/plot_top_virus_coverage.R \
        --vif_report {params.output_dir}{params.sample}.vif.init.filtered.abridged.tsv  \
        --virus_fai {params.virus_fai} \
        --bam {params.output_dir}{params.sample}.VirusDetect.sorted.igvjs.bam \
        --output_prefix {params.output_dir}{params.sample}.VirusDetect &&
        # make bed for igvjs
        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/create_igvjs_virus_bed.py \
        --summary {params.output_dir}{params.sample}.VirusDetect.virus_read_counts_summary.tsv \
        --output_prefix {params.output_dir}{params.sample}.VirusDetect \
        --num_top_viruses 20 &&
        python3 {params.dir_for_scripts}/nf-viral-starsolo/modified_files_from_viralintegration-0.1.1/bin/create_insertion_site_inspector_js.py \
        --VIF_summary_tsv {params.output_dir}{params.sample}.VirusDetect.igvjs.table.tsv \
        --json_outfile {params.output_dir}{params.sample}.VirusDetect.igvjs.json &&

        {params.dir_for_scripts}/nf-viral-starsolo/bamsifter/bamsifter \
        -c 100 \
        -o {params.output_dir}{params.sample}.VirusDetect.igvjs.reads.bam \
        {params.output_dir}{params.sample}.VirusDetect.sorted.igvjs.bam &&

        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/create_igvjs_virus_fa.py \
        {params.output_dir}{params.sample}.VirusDetect.igvjs.bed \
        {params.fasta_HPV}  \
        {params.output_dir}{params.sample}.VirusDetect.igvjs.fa &&

        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/make_VIF_igvjs_html.py \
        --html_template {params.dir_for_scripts}/nf-viral-starsolo/igvjs_VIF.html \
        --fusions_json {params.output_dir}{params.sample}.VirusDetect.igvjs.json \
        --input_file_prefix {params.output_dir}{params.sample}.VirusDetect.igvjs \
        --html_output {params.output_dir}{params.sample}.VirusDetect.igvjs.html &&

        python3 {params.dir_for_scripts}/nf-viral-starsolo/modified_files_from_viralintegration-0.1.1/bin/find_closest.py \
        -i {params.output_dir}{params.sample}.vif.refined.tsv \
        -o {params.output_dir}{params.sample}.vif.refined.wRefGeneAnnots.tsv \
        --gtf {params.gtf_file} &&
        python3 {params.dir_for_scripts}/nf-viral-starsolo/modified_files_from_viralintegration-0.1.1/bin/create_insertion_site_inspector_js.py \
            --VIF_summary_tsv {params.output_dir}{params.sample}.vif.refined.wRefGeneAnnots.tsv \
            --json_outfile {params.output_dir}{params.sample}.vif.igv.json &&

        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/region_gtf_to_bed.py \
            {params.output_dir}{params.sample}.vif.extract.gtf \
            > {params.output_dir}{params.sample}.vif.bed &&

        {params.dir_for_scripts}/nf-viral-starsolo/bamsifter/bamsifter \
            -c 100 \
            -o {params.output_dir}{params.sample}.vif.reads.bam \
            {params.output_dir}{params.sample}.vif.evidence.bam &&

        python3 {params.dir_for_scripts}/nf-viral-starsolo/viralintegration/bin/make_VIF_igvjs_html.py \
            --html_template {params.dir_for_scripts}/nf-viral-starsolo/igvjs_VIF.html \
            --fusions_json {params.output_dir}{params.sample}.vif.igv.json \
            --input_file_prefix {params.output_dir}{params.sample}.vif \
            --html_output {params.output_dir}{params.sample}.vif.html
        
        """

