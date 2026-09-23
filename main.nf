nextflow.enable.dsl = 2

include { MultiQC as MultiQC_Raw      } from './modules/multiqc.nf'
include { MultiQC as MultiQC_Trimmed  } from './modules/multiqc.nf'
include { MultiQC as MultiQC_Cleaned  } from './modules/multiqc.nf'

include { FastQC as FastQC_Raw         } from './modules/fastqc.nf'
include { FastQC as FastQC_Cleaned_R1  } from './modules/fastqc.nf'
include { FastQC as FastQC_Cleaned_R2  } from './modules/fastqc.nf'

include { FeatureCounts as FeatureCounts_unstranded          } from './modules/featurecounts.nf'
include { FeatureCounts as FeatureCounts_stranded            } from './modules/featurecounts.nf'
include { FeatureCounts as FeatureCounts_reversely_stranded  } from './modules/featurecounts.nf'



// =============================
// WORKFLOW
// =============================
workflow {

    // ------------------------------------------
    // STEP 1: CONCATENATE RAW FASTQ FILES
    // ------------------------------------------

    // Create a channel containing all FASTQ files that match the given glob pattern

     raw_reads_ch = channel
    // Match all relevant FASTQ extensions
    .fromPath("${params.raw_fastq_input_dir}/*.{fastq,fastq.gz,fq,fq.gz}")

    .map { file ->
        // Extract base name (without extension)
        def name = file.getBaseName().replaceFirst(/\.f(ast)?q(\.gz)?$/, '')
        // Remove .fastq, .fastq.gz, .fq, .fq.gz

        def tokens = name.split('_')
        def sample = tokens[0]                               // first token = sample name
        def read = (tokens.find { token -> token ==~ /(?i)^R[12]$/ } ?: 'R?').toUpperCase() // detect R1/R2 case-insensitively


        tuple(sample, "${sample}_${read}", file)
    }

    // Group by sample and read pair
    .groupTuple(by: [0, 1])

    // Optional: print what was grouped (useful for debugging)
    // .view { grouped_samples -> "Grouped samples: ${grouped_samples}" }

        
    // Concatenate files for each sample+read group
    // If a sample has multiple files (e.g., from different lanes), they will be concatenated to 1 file per R1/R2
    concat_fastq_out_ch = Concatenate_fastq(raw_reads_ch)
    
    // ------------------------------------------
    // STEP 2: FASTQC ON CONCATENATED FILES
    // ------------------------------------------

    raw_fastqc_out_ch = FastQC_Raw(concat_fastq_out_ch, 'fastqc_raw')

    // ------------------------------------------
    // STEP 3: MULTIQC ON RAW FASTQ FILES
    // ------------------------------------------

    raw_multiqc_input_ch = raw_fastqc_out_ch.fastqc_htmls // Collect all HTML reports
    .mix(raw_fastqc_out_ch.fastqc_zips)     // Mix with ZIP reports for MultiQC
    .collect() // Collect into a single channel for MultiQC
    .map { files -> tuple('multiqc_raw', files) } // Tag the collected files as "multiqc_raw"

    // Optional debug output
    // .view { raw_files -> "Raw MultiQC input files: ${raw_files}" }

    MultiQC_Raw(raw_multiqc_input_ch)
   

    // =======================================================
    // STEP 4: QC ON CONCATENATED FASTQ FILES WITH TRIM GALORE
    // =======================================================

    trim_galore_input_ch = concat_fastq_out_ch
    .groupTuple(by: 0, size: 2)  // (sample, [read_ids], [files])

    .map { sample, _reads, files -> // '_' before reads just means it's not being used
        // Match files dynamically by name (case-insensitive)
        def r1_file = files.find { file -> file.name =~ /(?i)_R1/ }
        def r2_file = files.find { file -> file.name =~ /(?i)_R2/ }

        tuple(sample, r1_file, r2_file)
    }

    // Keep only valid pairs
    .filter { _sample, r1_file, r2_file -> r1_file && r2_file } // '_' before sample just means it's not being used

    // Optional debug output
    // .view { item -> "Prepared for Trim Galore: ${item[0]} | R1: ${item[1].name} | R2: ${item[2].name}" }

    trim_galore_output_ch = TrimGalore(trim_galore_input_ch)

    // ------------------------------------------
    // STEP 5: MULTIQC ON TRIMMED FILES
    // ------------------------------------------

   trimmed_fastqc_ch = trim_galore_output_ch.fastqc_htmls
    // Get the FastQC HTML outputs from Trim Galore
        .mix(trim_galore_output_ch.fastqc_zips)
        // Combine them with the FastQC ZIP outputs
        .mix(trim_galore_output_ch.trimming_reports)
        // Mix in the trimming reports for MultiQC
        // .mix(trim_galore_output_ch.trimmed_reads_R1)
        // // Include the trimmed reads R1         
        // .mix(trim_galore_output_ch.trimmed_reads_R2)
        // Include the trimmed reads R2
        .flatten()
        // Convert grouped file lists into a flat stream of individual files
    
        .collect()
        // Gather all individual files into a single list for MultiQC
        .map { files -> tuple('multiqc_trimmed', files) } // Tag the collected files as "multiqc_raw"

        // Optional debug output
        // .view { trimmed -> "Trimmed MultiQC input files: ${trimmed}" }
       

    MultiQC_Trimmed(trimmed_fastqc_ch)
    // Run MultiQC on the collected FastQC outputs

    
    // ------------------------------------------------
    // STEP 6: SORTMERNA - remove rRNA from fastq files
    // ------------------------------------------------

    // Use the emited tuple with trimmed reads from Trim Galore as input for SortMeRNA
    sortmerna_input_ch = trim_galore_output_ch.trimmed_reads

    // Optional debug output
    // sortmerna_input_ch.view { item -> "SortMeRNA input: ${item[0]} -> ${item[1].getName()}, ${item[2].getName()}" }

    sortmerna_output_ch = SortMeRNA(sortmerna_input_ch)

    // ------------------------------------------------
    // STEP 7: Split interleaved reads from SortMeRNA
    // ------------------------------------------------

    SplitCleanReads(sortmerna_output_ch.sortmerna_clean)

    // ------------------------------------------------
    // STEP 8: FastQC on splitted cleaned fastq files
    // ------------------------------------------------

    cleaned_fastqc_out_R1_ch = FastQC_Cleaned_R1(SplitCleanReads.out.trimmed_reads_R1, 'fastqc_cleaned')
    cleaned_fastqc_out_R2_ch = FastQC_Cleaned_R2(SplitCleanReads.out.trimmed_reads_R2, 'fastqc_cleaned')

    // ------------------------------------------
    // STEP 9: MULTIQC ON CLEANED FILES
    // ------------------------------------------

    cleaned_multiqc_input_ch = cleaned_fastqc_out_R1_ch.fastqc_htmls
        .mix(cleaned_fastqc_out_R2_ch.fastqc_htmls)
        .mix(cleaned_fastqc_out_R1_ch.fastqc_zips)
        .mix(cleaned_fastqc_out_R2_ch.fastqc_zips)
        .collect()
        .map { files -> tuple('multiqc_cleaned', files) }

    MultiQC_Cleaned(cleaned_multiqc_input_ch)


    // ------------------------------------------
    // STEP 10: STAR ALIGNMENT
    // ------------------------------------------

    // Create channels for the reference genome and GTF files
    reference_fasta_ch = channel.fromPath(params.reference_genome)
    reference_gtf_ch   = channel.fromPath(params.reference_gtf)

    // Build the index (runs once) and get its emitted channel
    genome_idx_ch = GenomeGenerate(reference_fasta_ch, reference_gtf_ch)

    // Join the R1 and R2 files on their sample prefix
    align_reads_input_ch = SplitCleanReads.out.trimmed_reads_R1          // channel #1
    .join(SplitCleanReads.out.trimmed_reads_R2)                          // channel #2
    .combine(genome_idx_ch)                                              // add the genome index channel to each tuple
   
    .map { id, _r1, pathR1, _r2, pathR2, idx_dir ->                          // destructure, prefix '_' in r1 and r2 means they are not being used
        tuple(id, pathR1, pathR2, idx_dir)                                   // (sample, R1, R2)
    }

    // Optional debug output
    // .view { item -> "align reads input ch: ${item}" } 
       
        
    // Input cleaned fastq files from sortmerna
    alignreads_output_ch = AlignReads(align_reads_input_ch)

    // Optional debug output
    // alignreads_output_ch.view { item -> "alignreads_output_ch: ${item}"}

    // ------------------------------------------
    // STEP 11: MarkDuplicates
    // ------------------------------------------

    markduplicates_output_ch = MarkDuplicates(alignreads_output_ch)
    .map { sample_id, bam_file, _bai_file -> // prefix '_' before bai_file suppress warning that it's not being used
        tuple(sample_id, bam_file) 
    }

    
    // ------------------------------------------
    // STEP 12: SplitNCigarReads
    // ------------------------------------------

    splitncigarreads_output_ch = SplitNCigarReads(markduplicates_output_ch)

    // Optional debug output
    // splitncigarreads_output_ch.view { item -> "splitncigarreads_output_ch: ${item}"}

    // ------------------------------------------
    // STEP 13: BaseRecalibrator
    // ------------------------------------------

    base_recalibrator_out_ch = BaseRecalibrator(splitncigarreads_output_ch)

    // Optional debug output
    // base_recalibrator_out_ch.view { item -> "base_recalibrator_out_ch: ${item}"}

    // Remap to pass BAM + recal table to ApplyBQSR
    apply_bqsr_input_ch = base_recalibrator_out_ch
        .join(splitncigarreads_output_ch) // .join is best practice for joining channels when they are related
        .map { sample_id, recal_table, bam_file ->
            tuple(sample_id, bam_file, recal_table)
        }

        // Optional debug output
        // .view { item -> "apply_bqsr_input_ch: ${item}" }


    ApplyBQSR(apply_bqsr_input_ch)

    // ------------------------------------------
    // STEP 14: Haplotypecaller
    // ------------------------------------------

    HaplotypeCaller(ApplyBQSR.out)

    // ------------------------------------------
    // STEP 15: VariantFiltration
    // ------------------------------------------

    VariantFiltration(HaplotypeCaller.out)

    // ------------------------------------------
    // STEP 16: WASP step 1: ExtractVcfSNPs
    // ------------------------------------------
    
    ExtractVcfSNPs(VariantFiltration.out)

    // ------------------------------------------
    // STEP 17: WASP step 2: FindIntersectingSNPs
    // ------------------------------------------

    findintersectingsnps_input_ch = ApplyBQSR.out
    .join(ExtractVcfSNPs.out)
    // Optional debug output
    // .view {before_map -> "Before map: ${before_map}"}
    .map {sample_id_bam, bam_file, txt_snps ->
    tuple(sample_id_bam, bam_file, txt_snps)
    }
    // Optional debug output
    // .view {after_map -> "findintersectingsnps_input_ch: ${after_map}"}

    FindIntersectingSNPs(findintersectingsnps_input_ch)

    // ------------------------------------------
    // STEP 18: WASP step 3: RemapReads
    // ------------------------------------------

    remap_reads_input_ch = FindIntersectingSNPs.out
        .map { sample_id, fq1, fq2, _single, _keep_bam, _to_remap_bam -> // prefix '_' used to supress warning that those items are not being used 
        tuple(sample_id, fq1, fq2) }
        // Optional debug output
        // .view { item -> "remap_reads_input: ${item}"}

    RemapReads(remap_reads_input_ch)

    // ------------------------------------------
    // STEP 19: WASP step 4: FilterRemapReads
    // ------------------------------------------

    filter_remap_reads_input_ch = FindIntersectingSNPs.out
        .join(RemapReads.out)
        // Optional debug output
        // .view {before_map -> "Before map: ${before_map}"}
        .map { sample_id, _fq1, _fq2, _single, _keep_bam, to_remap_bam, remap_bam -> // prefix '_' used to supress warning that those items are not being used
        tuple(sample_id, to_remap_bam, remap_bam) }
        // Optional debug output
        // .view {after_map -> "After map: ${after_map}"}

        FilterRemapReads(filter_remap_reads_input_ch)

    // ------------------------------------------
    // STEP 20: WASP step 4: MergeMappedReads
    // ------------------------------------------

    merge_remap_reads_input_ch = FindIntersectingSNPs.out
        .join(FilterRemapReads.out)
        // Optional debug output
        // .view { item -> "Before map: ${item}"}
        .map { sample_id, _fq1, _fq2, _single, keep_bam, _to_remap_bam, remapped_filtered_bam -> // prefix '_' used to supress warning that those items are not being used
        tuple(sample_id, keep_bam, remapped_filtered_bam) }
        // Optional debug output
        // .view {item -> "After map: ${item}"}


    MergeRemapReads(merge_remap_reads_input_ch)

    // ------------------------------------------
    // STEP 21: WASP step 5: FilterDuplicateReads 
    // ------------------------------------------

    FilterDuplicateReads(MergeRemapReads.out)

    // ------------------------------------------
    // STEP 22: WASP step 6: SelectBiallelicSites 
    // ------------------------------------------

    SelectBiallelicSites(VariantFiltration.out)

    // ------------------------------------------
    // STEP 23: WASP step 7: ASEReadCounter 
    // ------------------------------------------

    ase_read_counter_input_ch = FilterDuplicateReads.out
        .join(SelectBiallelicSites.out)
        // Optional debug output
        // .view { item -> "before map: ${item}"}
        .map { sample_id, dedup_sort_bam, _dedup_sort_bai, biallelic_variants -> // prefix '_' used to supress warning that those items are not being used
        tuple(sample_id, dedup_sort_bam, biallelic_variants) }
        // Optional debug output
        // .view { item -> "after map: ${item}"}

    ASEReadCounter(ase_read_counter_input_ch)

    // ------------------------------------------
    // STEP 24: FeatureCounts 
    // ------------------------------------------

    // Perform featureCounts on the deduplicated BAM files from FilterDuplicateReads, same as used for ASEReadCounter
    featurecounts_input_ch = FilterDuplicateReads.out
        .map {_sample_id, dedup_sort_bam, _dedup_sort_bai -> dedup_sort_bam } // prefix '_' used to supress warning that those items are not being used
        .collect()


    FeatureCounts_unstranded(featurecounts_input_ch, 0)
    FeatureCounts_stranded(featurecounts_input_ch, 1)
    FeatureCounts_reversely_stranded(featurecounts_input_ch, 2)

}



// =============================
// PROCESS: CONCATENATE FASTQ
// =============================
process Concatenate_fastq {

    tag "${sample_read}"

    input:
    tuple val(sample), val(sample_read), path(fastq_files)

    output:
    tuple val(sample), val(sample_read), path("${sample_read}_${task.process}.fastq.gz")
    

    script:
    // This is groovy code and therefore outside the triple quotes
    // Define the output file name based on the sample name
    def output_file = "${sample_read}_${task.process}.fastq.gz"
       
    // - If there's only one input FASTQ file, just copy it (no need to concatenate).
    // - If there are multiple FASTQ files (e.g., multiple lanes for the same sample),
    //   concatenate them into a single output file using 'cat'.
    def cmd = fastq_files.size() == 1
        ? "cp ${fastq_files[0]} ${output_file}"                        // Case: only one file — use 'cp'
        : "cat ${fastq_files.join(' ')} > ${output_file}"              // Case: multiple files — use 'cat' to merge

    

    """

    # Enable strict Bash settings for safety:
    # -e  : exit immediately if any command fails
    # -u  : treat unset variables as errors
    # -o pipefail : fail if any command in a pipeline fails (not just the last one)
    set -euo pipefail

    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] ${task.process} lanes for: ${task.tag}" >> "\$PWD/.command.out"
    
    if ! ${cmd}; then
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
}


// =============================
// PROCESS: TrimGalore
// =============================

process TrimGalore {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple  val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads

    path   "*_fastqc.html",            emit: fastqc_htmls
    path   "*_fastqc.zip",             emit: fastqc_zips
    path   "*trimming_report.txt",     optional: true, emit: trimming_reports

    script:
     
    """

    # Enable strict Bash settings for safety:
    # -e  : exit immediately if any command fails
    # -u  : treat unset variables as errors
    # -o pipefail : fail if any command in a pipeline fails (not just the last one)
    set -euo pipefail

    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"
   
    echo "🔍 Running ${task.process} on: ${read1} and ${read2} for sample ${task.tag}"

    # Optional trimming parameters: -a " AGATCGGAAGAGC -a G{50}" -a2 " AGATCGGAAGAGC -a G{50}"

    if ! trim_galore \
        --paired \
        --quality 20 \
        --length 20 \
        --cores 4 \
        --gzip \
        --fastqc \
        --output_dir . \
        "${read1}" "${read2}" \
        1>> \$PWD/.command.out 2>&1; then
        echo "❌ ${task.process} failed for ${read1} and ${read2} for sample ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} was succesful for ${read1} and ${read2} for sample ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
}


// =============================
// PROCESS: SortMeRNA
// =============================

process SortMeRNA {

    publishDir "${params.sortmerna_output_dir}/${task.tag}", mode: params.publish_mode

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_merged_clean.fq.gz"), emit: sortmerna_clean
    // path("${sample_id}_rRNA.fq.gz"), emit: sortmerna_rrna
    
    script:

    """

    # Enable strict Bash settings for safety:
    # -e  : exit immediately if any command fails
    # -u  : treat unset variables as errors
    # -o pipefail : fail if any command in a pipeline fails (not just the last one)
    set -euo pipefail

    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "🔍 Running ${task.process} on: ${read1} and ${read2} for sample ${task.tag}"

    if ! sortmerna \
        --ref ${params.sortmerna_ref} \
        --idx-dir ${params.sortmerna_index_dir} \
        --reads ${read1} \
        --reads ${read2} \
        --paired_in --fastx \
        --other ${sample_id}_merged_clean \
        --aligned ${sample_id}_rRNA \
        --threads ${task.cpus} \
        --workdir ${sample_id}_tmp \
        1>> \$PWD/.command.out 2>&1; then
        echo "❌ ${task.process} failed for ${read1} and ${read2} for sample ${task.tag}" >> \$PWD/.command.err
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} was successful for ${read1} and ${read2} for sample ${task.tag}" >> \$PWD/.command.out
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
}


// ==========================
// PROCESS: SplitCleanReads 
// ==========================

process SplitCleanReads {

    publishDir "${params.sortmerna_output_dir}/${task.tag}", mode: params.publish_mode

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(clean_reads)

    output:
    tuple val(sample_id), val("${sample_id}_R1"), path("${sample_id}_clean_R1.fastq.gz"), emit: trimmed_reads_R1
    tuple val(sample_id), val("${sample_id}_R2"), path("${sample_id}_clean_R2.fastq.gz"), emit: trimmed_reads_R2

    // tuple val(sample_id), path("${sample_id}_clean_R1.fastq.gz"), path("${sample_id}_clean_R2.fastq.gz")

    script:

    """

    # Enable strict Bash settings for safety:
    # -e  : exit immediately if any command fails
    # -u  : treat unset variables as errors
    # -o pipefail : fail if any command in a pipeline fails (not just the last one)
    set -euo pipefail

    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "🔍 Running ${task.process} to split interleaved reads for sample ${task.tag}"

    if ! reformat.sh \
        in=${clean_reads} \
        out1=${sample_id}_clean_R1.fastq.gz \
        out2=${sample_id}_clean_R2.fastq.gz \
        1>> \$PWD/.command.out 2>&1; then
        echo "❌ ${task.process} failed for ${clean_reads} for sample ${task.tag}" >> \$PWD/.command.err
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} was successful for ${clean_reads} for sample ${task.tag}" >> \$PWD/.command.out
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
}



// ============================
// PROCESS: STAR genomeGenerate 
// ============================

process GenomeGenerate {

    publishDir params.star_genome, mode: params.publish_mode   // optional: copy final index to your desired location

    tag "STAR genome index"

    input:
    path reference_fasta
    path reference_gtf

    output:
    path "genome_index", emit: star_index_dir     // emit the index dir produced in workdir

    script:
    """

    # Enable strict Bash settings for safety:
    # -e  : exit immediately if any command fails
    # -u  : treat unset variables as errors
    # -o pipefail : fail if any command in a pipeline fails (not just the last one)
    set -euo pipefail


    # Create the genome index directory if it doesn't exist inside the work directory
    mkdir -p genome_index

    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process}" >> "\$PWD/.command.out"

    if ! STAR \
        --runMode genomeGenerate \
        --genomeDir genome_index \
        --genomeFastaFiles ${reference_fasta} \
        --sjdbGTFfile ${reference_gtf} \
        --sjdbOverhang 100 \
        --runThreadN ${task.cpus} \
        1>> \$PWD/.command.out 2>&1; then
        echo "❌ ${task.process} failed" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh" "*_Log.out" "*_Log.final.out"
        exit 1
    else
        echo "✅ ${task.process} done" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh" "*_Log.out" "*_Log.final.out"
        exit 0
    fi
    """
}


// ==========================
// PROCESS: STAR alignReads 
// ==========================


process AlignReads {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(read1), path(read2), path(star_index_dir)
    

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    

    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # ------------------------------------------------------------------
    # 1. Derive read-group (RG) fields from the first FASTQ header
    #    Header template: @<instrument>:<run>:<FLOWCELL>:<LANE>:...
    # ------------------------------------------------------------------
    
    header=\$(zcat -f ${read1} | head -1)                               # @NB501805:7:H3GG7BGX3:1:11101:...
    flowcell=\$(echo "\$header" | cut -d ':' -f3)                       # e.g. H3GG7BGX3
    lane=\$(echo     "\$header" | cut -d ':' -f4)                       # e.g. 1
    barcode=\$(echo  "\$header" | awk '{print \$2}' | cut -d ':' -f4)   # e.g. CGATGT

    LB="\${flowcell}_\${barcode}"       # library
    PU="\${flowcell}.\${lane}"          # platform-unit
    PL="ILLUMINA"                       # platform
    RG_LINE="ID:${sample_id} SM:${sample_id} LB:\${LB} PU:\${PU} PL:\${PL}"

    echo "[DEBUG] RG → \$RG_LINE" >> "\$PWD/.command.out"

    # Run STAR and catch failures manually
    if ! STAR \
        --runMode alignReads \
        --genomeDir ${star_index_dir} \
        --runThreadN ${task.cpus} \
        --readFilesIn ${read1} ${read2} \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --sjdbOverhang 100 \
        --outFilterMismatchNoverLmax 0.3 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample_id}_ \
        --outSAMattrRGline \$RG_LINE \
        1>> \$PWD/.command.out 2>&1; then
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh" "*_Log.out" "*_Log.final.out"
        exit 1
    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh" "*_Log.out" "*_Log.final.out"
        exit 0
    fi
    """
}

// ============================
// PROCESS: STAR alignReadsWASP
// ============================


process AlignReadsWASP {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    

    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # ------------------------------------------------------------------
    # 1. Derive read-group (RG) fields from the first FASTQ header
    #    Header template: @<instrument>:<run>:<FLOWCELL>:<LANE>:...
    # ------------------------------------------------------------------
    
    header=\$(zcat -f ${read1} | head -1)                               # @NB501805:7:H3GG7BGX3:1:11101:...
    flowcell=\$(echo "\$header" | cut -d ':' -f3)                       # e.g. H3GG7BGX3
    lane=\$(echo     "\$header" | cut -d ':' -f4)                       # e.g. 1
    barcode=\$(echo  "\$header" | awk '{print \$2}' | cut -d ':' -f4)   # e.g. CGATGT

    LB="\${flowcell}_\${barcode}"       # library
    PU="\${flowcell}.\${lane}"          # platform-unit
    PL="ILLUMINA"                       # platform
    RG_LINE="ID:${sample_id} SM:${sample_id} LB:\${LB} PU:\${PU} PL:\${PL}"

    echo "[DEBUG] RG → \$RG_LINE" >> "\$PWD/.command.out"

    # Run STAR and catch failures manually
    if ! STAR \
        --runMode alignReads \
        --genomeDir ${params.star_genome} \
        --genomeLoad NoSharedMemory \
        --runThreadN ${task.cpus} \
        --readFilesIn ${read1} ${read2} \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --alignEndsType EndToEnd \
        --outFilterMultimapNmax 1 \
        --sjdbOverhang 100 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --outFilterMismatchNoverLmax 0.04 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM NM MD jM jI rB MC vA vG vW \
        --waspOutputMode SAMtag \
        --outFileNamePrefix ${sample_id}_ \
        --outSAMattrRGline \$RG_LINE \
        --outSAMunmapped Within \
        1>> \$PWD/.command.out 2>&1; then
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh" "*_Log.out" "*_Log.final.out"
        exit 1
    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh" "*_Log.out" "*_Log.final.out"
        exit 0
    fi
    """
}
    
// ==========================
// PROCESS: MarkDuplicates
// ==========================

process MarkDuplicates {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("${sample_id}_MarkDuplicates.bam"), path("${sample_id}_MarkDuplicates.bai")

    script:
    

    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Run Picard MarkDuplicates
    if ! java -Xmx20g -jar /ludc/Home/jonas_a/.conda/pkgs/picard-3.3.0-hdfd78af_0/share/picard-3.3.0-0/picard.jar MarkDuplicates \
        I=${input_bam} \
        O=${sample_id}_MarkDuplicates.bam \
        M=${sample_id}_marked-dup-metrics.txt \
        CREATE_INDEX=true \
        ASSUME_SORTED=true \
        TMP_DIR=${params.tmp_dir} \
        VALIDATION_STRINGENCY=SILENT \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
    
}


// ==========================
// PROCESS: SplitNCigarReads
// ==========================

process SplitNCigarReads {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("${sample_id}_SplitNCigarReads.bam")

    script:
    

    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Run GATK SplitNCigarReads
    if ! gatk SplitNCigarReads \
        -R ${params.reference_genome} \
        -I ${input_bam} \
        -O ${sample_id}_SplitNCigarReads.bam \
        1>> \$PWD/.command.out 2>&1; then
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
    
}

// ==========================
// PROCESS: BaseRecalibrator
// ==========================

process BaseRecalibrator {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("${sample_id}_recal_data.table")

    script:
    

    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Run GATK BaseRecalibrator
    if ! java -Djava.io.tmpdir=${params.tmp_dir} \
        -XX:GCTimeLimit=50 \
        -XX:GCHeapFreeLimit=10 \
        -XX:+PrintFlagsFinal \
        -Xlog:gc*:file=gc_log.log:time,uptime,level,tags \
        -Xms4000m \
        -jar /ludc/Home/jonas_a/.conda/envs/variant_calling_py310/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
        BaseRecalibrator \
        -R ${params.reference_genome} \
        -I ${input_bam} \
        -OQ \
        --known-sites ${params.snp_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
        --known-sites ${params.snp_sites}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        --known-sites ${params.snp_sites}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -O "${sample_id}_recal_data.table" \
        1>> \$PWD/.command.out 2>&1; then
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
    
}

// ================================
// PROCESS: Apply BaseRecalibration
// ================================

process ApplyBQSR {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(input_bam), path(recal_table)

    output:
    tuple val(sample_id), path("${sample_id}_recalibrated.bam")

    script:
    

    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Run GATK BaseRecalibrator
    if ! java -Djava.io.tmpdir=${params.tmp_dir} \
    -XX:GCTimeLimit=50 \
    -XX:GCHeapFreeLimit=10 \
    -XX:+PrintFlagsFinal \
    -Xlog:gc*:file=gc_log.log:time,uptime,level,tags \
    -Xms4000m \
    -jar /ludc/Home/jonas_a/.conda/envs/variant_calling_py310/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
    ApplyBQSR \
    --add-output-sam-program-record \
    --use-original-qualities \
    --bqsr-recal-file ${recal_table} \
    -R ${params.reference_genome} \
    -I ${input_bam} \
    -O "${sample_id}_recalibrated.bam" \
    1>> \$PWD/.command.out 2>&1; then

        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
    
}


// ========================
// PROCESS: HaplotypeCaller
// ========================

process HaplotypeCaller {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("${sample_id}_haplotypecaller.vcf.gz")

    script:
    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Run GATK HaplotypeCaller
    if ! java -Djava.io.tmpdir=${params.tmp_dir} \
        -Xms6000m \
        -XX:GCTimeLimit=50 \
        -XX:GCHeapFreeLimit=10 \
        -XX:+PrintFlagsFinal \
        -Xlog:gc*:file=gc_log.log:time,uptime,level,tags \
        -jar /ludc/Home/jonas_a/.conda/envs/variant_calling_py310/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
        HaplotypeCaller \
        -R ${params.reference_genome} \
        -I ${input_bam} \
        -O "${sample_id}_haplotypecaller.vcf.gz" \
        --dont-use-soft-clipped-bases \
        --dbsnp ${params.snp_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1

    else       
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
}

// ==========================
// PROCESS: VariantFiltration
// ==========================

process VariantFiltration {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(input_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz")

    script:
    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Index the input VCF
    if ! java -Djava.io.tmpdir=${params.tmp_dir} \
        -Xms6000m \
        -XX:GCTimeLimit=50 \
        -XX:GCHeapFreeLimit=10 \
        -XX:+PrintFlagsFinal \
        -Xlog:gc*:file=gc_log.log:time,uptime,level,tags \
        -jar /ludc/Home/jonas_a/.conda/envs/variant_calling_py310/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
        IndexFeatureFile -I ${input_vcf} \
        1>> \$PWD/.command.out 2>&1; then

            echo "❌ IndexFeatureFile failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ IndexFeatureFile done for ${task.tag}" >> "\$PWD/.command.out"
    fi

    # Run GATK VariantFiltration
    if ! java -Djava.io.tmpdir=${params.tmp_dir} \
        -Xms6000m \
        -XX:GCTimeLimit=50 \
        -XX:GCHeapFreeLimit=10 \
        -XX:+PrintFlagsFinal \
        -Xlog:gc*:file=gc_log.log:time,uptime,level,tags \
        -jar /ludc/Home/jonas_a/.conda/envs/variant_calling_py310/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
        VariantFiltration \
        -R ${params.reference_genome} \
        -V ${input_vcf} \
        -window 35 \
        -cluster 3 \
        --filter-name "FS" -filter "FS > 30.0" \
        --filter-name "QD" -filter "QD < 2.0" \
        -O "${sample_id}_filtered.vcf.gz" \
        1>> \$PWD/.command.out 2>&1; then
        
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
}

// =============================
// PROCESS: ExtractVcfSNPs
// =============================

process ExtractVcfSNPs {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(input_vcf)

    output:
    tuple val(sample_id), path("${sample_id}") // folder with per-contig SNP txt.gz files

    script:
    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Create output folder
    mkdir -p ${sample_id}

    # Extract SNPs and split by contig (column 1)
    if ! zcat ${input_vcf} | awk '!/^#/ {print > "'${sample_id}'/"\$1".tmp"}'; then
        echo "❌ Failed splitting SNPs for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    fi

    # Format and compress: POS, REF, ALT — no header
    for file in ${sample_id}/*.tmp; do
        contig=\$(basename "\$file" .tmp)
        awk '{print \$2, \$4, \$5}' "\$file" | gzip > "${sample_id}/\${contig}.snps.txt.gz"
        rm "\$file"
    done

    echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
    move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
    exit 0
    """
}

// =============================
// PROCESS: FindIntersectingSNPs
// =============================

process FindIntersectingSNPs {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(input_bam), path(txt_snps)

    output:
    tuple val(sample_id), path("${sample_id}_recalibrated.remap.fq1.gz"), path("${sample_id}_recalibrated.remap.fq2.gz"), path("${sample_id}_recalibrated.remap.single.fq.gz"), path("${sample_id}_recalibrated.keep.bam"), path("${sample_id}_recalibrated.to.remap.bam")

    script:
    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Run find_intersecting_snps.py
    if ! python ${params.wasp_dir}/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir . \
        --snp_dir ${txt_snps} \
        ${input_bam} \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1

    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
}

// =============================
// PROCESS: RemapReads
// =============================

process RemapReads {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}_remapped_Aligned.sortedByCoord.out.bam")

    script:
    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔁 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Run STAR for remapping
    if ! STAR \
        --runMode alignReads \
        --genomeDir "${params.star_genome}/genome_index" \
        --runThreadN ${task.cpus} \
        --readFilesIn ${fq1} ${fq2} \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --sjdbOverhang 100 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "${sample_id}_remapped_" \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1

    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
}

// =============================
// PROCESS: FilterRemapReads
// =============================

process FilterRemapReads {

    tag { "${sample_id}" }

    // recalibrated_to_remap: input BAM file containing original set of reads that needed to be remapped after having their alleles  flipped. This file is output by the find_intersecting_snps.py script.
    // remap_bam: input BAM file containing remapped reads (with flipped alleles)
    input:
    tuple val(sample_id), path(recalibrated_to_remap), path(remapped_bam)

    output:
    tuple val(sample_id), path("${sample_id}_remapped_filtered_keep.bam")

    script:
    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔁 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Run STAR for remapping
    if ! python ${params.wasp_dir}/mapping/filter_remapped_reads.py \
        ${recalibrated_to_remap} \
        ${remapped_bam} \
		${sample_id}_remapped_filtered_keep.bam \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1

    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
}

// =============================
// PROCESS: MergeRemappedReads
// =============================

process MergeRemapReads {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(recalibrated_to_keep), path(remapped_filtered_bam_to_keep)

    output:
    tuple val(sample_id), path("${sample_id}_keep_merged_sorted.bam"), path("${sample_id}_keep_merged_sorted.bam.bai")

    script:
    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔁 Merging bam for sample: ${task.tag}" >> "\$PWD/.command.out"


    # Step 1: Merge BAMs
    if ! samtools merge \
        ${sample_id}_keep_merged.bam \
        ${recalibrated_to_keep} \
        ${remapped_filtered_bam_to_keep} \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ Merge failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1

    else
        echo "✅ Merge done for ${task.tag}" >> "\$PWD/.command.out"   
    fi

    echo "[DEBUG] 🔁 Sorting bam for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Step 2: Sort merged BAM
    if ! samtools sort \
    -o ${sample_id}_keep_merged_sorted.bam \
       ${sample_id}_keep_merged.bam \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ Sort failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    
    else
        echo "✅ Sort done for ${task.tag}" >> "\$PWD/.command.out"
    fi

    echo "[DEBUG] 🔁 Indexing  bam for sample: ${task.tag}" >> "\$PWD/.command.out"
    
    # Step 3: Index sorted BAM
    if ! samtools index \
        ${sample_id}_keep_merged_sorted.bam \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ Indexing failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    
    else
        echo "✅ Indexing done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
    fi
    """
}

// =============================
// PROCESS: FilterDuplicateReads
// =============================

process FilterDuplicateReads {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(keep_merged_sorted_bam), path(keep_merged_sorted_bam_bai)

    output:
    tuple val(sample_id), path("${sample_id}_dedup_sort.bam"), path("${sample_id}_dedup_sort.bam.bai")

    script:
    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔁 Merging bam for sample: ${task.tag}" >> "\$PWD/.command.out"


    # Step 1: Remove duplicate reads
    if ! python ${params.wasp_dir}/mapping/rmdup_pe.py \
        ${keep_merged_sorted_bam} \
        ${sample_id}_dedup.bam \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ Remove duplicate reads failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1

    else
        echo "✅ Remove duplicate reads done for ${task.tag}" >> "\$PWD/.command.out"   
    fi

    echo "[DEBUG] 🔁 Sorting bam for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Step 2: Sort deduped BAM
    if ! samtools sort \
    -o ${sample_id}_dedup_sort.bam \
       ${sample_id}_dedup.bam \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ Sort failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    
    else
        echo "✅ Sort done for ${task.tag}" >> "\$PWD/.command.out"
    fi

    echo "[DEBUG] 🔁 Indexing  bam for sample: ${task.tag}" >> "\$PWD/.command.out"
    
    # Step 3: Index sorted BAM
    if ! samtools index \
        ${sample_id}_dedup_sort.bam \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ Indexing failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    
    else
        echo "✅ Indexing done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
    fi
    """
}

// =============================
// PROCESS: SelectBiallelicSites
// =============================

process SelectBiallelicSites {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(filtered_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_biallelic_variants.vcf.gz")

    script:
    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔁 ${task.process} bam for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Index the input VCF
    if ! java -Djava.io.tmpdir=${params.tmp_dir} \
        -Xms6000m \
        -XX:GCTimeLimit=50 \
        -XX:GCHeapFreeLimit=10 \
        -XX:+PrintFlagsFinal \
        -Xlog:gc*:file=gc_log.log:time,uptime,level,tags \
        -jar /ludc/Home/jonas_a/.conda/envs/variant_calling_py310/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
        IndexFeatureFile -I ${filtered_vcf} \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ IndexFeatureFile failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ IndexFeatureFile done for ${task.tag}" >> "\$PWD/.command.out"
    fi

    # Select biallelic sites
    if ! java -Djava.io.tmpdir=${params.tmp_dir} \
        -jar /ludc/Home/jonas_a/.conda/envs/variant_calling_py310/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
        SelectVariants \
        -R ${params.reference_genome} \
        -V ${filtered_vcf}\
        --select-type-to-include SNP \
        --restrict-alleles-to BIALLELIC \
        -O "${sample_id}_biallelic_variants.vcf.gz" \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1

    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0 
    fi
    """
}

// =============================
// PROCESS: AseReadCounter
// =============================

process ASEReadCounter {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(wasp_sort_dedup_bam), path(wasp_biallelic_variants)

    output:
    tuple val(sample_id), path("${sample_id}_ASE_read_count.table")

    script:
    """
    # Source shared functions
    source "${params.modules_dir}/nextflow_functions.sh"

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    # Index the input VCF
    if ! java -Djava.io.tmpdir=${params.tmp_dir} \
        -Xms6000m \
        -XX:GCTimeLimit=50 \
        -XX:GCHeapFreeLimit=10 \
        -XX:+PrintFlagsFinal \
        -Xlog:gc*:file=gc_log.log:time,uptime,level,tags \
        -jar /ludc/Home/jonas_a/.conda/envs/variant_calling_py310/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
        IndexFeatureFile -I ${wasp_biallelic_variants} \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ IndexFeatureFile failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ IndexFeatureFile done for ${task.tag}" >> "\$PWD/.command.out"
    fi

    echo "[DEBUG] 🔁 ${task.process} bam for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Select biallelic sites
    if ! java -Djava.io.tmpdir=${params.tmp_dir} \
        -jar /ludc/Home/jonas_a/.conda/envs/variant_calling_py310/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
        ASEReadCounter \
        -R ${params.reference_genome} \
        -I ${wasp_sort_dedup_bam} \
        -V ${wasp_biallelic_variants} \
        -O ${sample_id}_ASE_read_count.table \
        1>> \$PWD/.command.out 2>&1; then

        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1

    else
        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0 
    fi
    """
}
