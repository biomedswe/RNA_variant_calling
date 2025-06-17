nextflow.enable.dsl = 2


// =============================
// WORKFLOW
// =============================
workflow {

    // ------------------------------------------
    // STEP 1: CONCATENATE RAW FASTQ FILES
    // ------------------------------------------
    Channel
        .fromPath("${params.raw_fastq_input_dir}/*.fastq.gz")  // Step 1: Read FASTQ files from input directory, include only files with .fastq.gz extension

        // Step 2: Validate filenames with expected pattern
        .filter { file ->
            def name = file.getName()
            return name ==~ /^HTL.*\.fastq\.gz$/ // Ensure the file name starts with "HTL" and ends with ".fastq.gz"
        }

        // Step 3: Read FASTQ files and extract a clean sample + read name
        .map { file ->
            def name = file.getBaseName()  // Removes the ".fastq.gz" suffix
            def tokens = name.split('_')  // Split on underscores
            def sample = tokens[0]        // First token = sample name
            def read = tokens.find { it ==~ /^R[12]$/ } ?: 'R?'  // Read pair
            tuple("${sample}_${read}", file) // Group by sample + read
        }

        // Step 4: Group files by sample+read (e.g. HTL123_R1)
        .groupTuple()

       

        | concat_fastq  // Step 5: Concatenate files for each sample+read group


    // ------------------------------------------
    // STEP 2: FASTQC ON CONCATENATED FILES
    // ------------------------------------------
    Channel
        .fromPath("${params.concat_fastq_dir}/*.fastq.gz")
        | fastqc

    // -----------------------------------------------
    // STEP 3: TRIM GALORE on concatenated FASTQ files
    // -----------------------------------------------   

    // Channel
    // .fromFilePairs("${params.trimgalore_input_dir}/*_R{1,2}.fastq.gz", flat: true)
    // .filter { sample_id, reads -> reads.size() == 2 }
    // .map { sample_id, reads ->
    //     println "[DEBUG] Matched sample: ${sample_id}"
    //     println "[DEBUG] → R1: ${reads[0].getName()}"
    //     println "[DEBUG] → R2: ${reads[1].getName()}"
    //     tuple(sample_id, reads[0], reads[1])
    // }
    // .set { trim_galore_input_ch }

    //  trim_galore(trim_galore_input_ch)

    Channel
        .fromFilePairs("${params.trimgalore_input_dir}/*_R{1,2}.fastq.gz", flat: false) // flat: false → files are symlinks to nextflow working directory and not in the same directory, so false is needed
        .filter { sample_id, reads -> reads.size() == 2 }
        .map    { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
    | trim_galore

    // ------------------------------------------
    // STEP 4: PICARD FastqToSam (paired-end)
    // ------------------------------------------
    // Channel
    //     .fromFilePairs("${params.trimgalore_output_dir}/*_R{1,2}.fastq.gz", flat: true)
    //     .filter { sample_id, reads -> reads.size() == 2 }
    //     .map    { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
    //     | fastqtosam
}



// =============================
// PROCESS: CONCATENATE FASTQ
// =============================
process concat_fastq {

    tag "${sample_read}"

    input:
    tuple val(sample_read), path(fastq_files)

    output:
    path("${sample_read}.fastq.gz")

    publishDir params.concat_fastq_dir, mode: params.publish_mode

    script:
    def output_file = "${sample_read}.fastq.gz"
    def cmd = fastq_files.size() == 1 
        ? "cp ${fastq_files[0]} ${output_file}" 
        : "cat ${fastq_files.join(' ')} > ${output_file}"

    """
    mkdir -p ${params.log_dir}
    echo "→ Merging files for ${sample_read}" >&2
    ${cmd}
    echo "✅ Created: ${output_file}" >&2
    """
}


// =============================
// PROCESS: FASTQC
// =============================
process fastqc {

    tag { input_file.getBaseName() }

    input:
    path input_file

    output:
    path "*.html"
    path "*.zip"

    publishDir params.fastqc_dir, mode: params.publish_mode

    script:
    """
    echo "🔍 Running FastQC on ${input_file.getName()}" >&2
    fastqc ${input_file} --outdir .
    echo "✅ FastQC done for ${input_file.getName()}" >&2
    """
}

// =============================
// PROCESS: TrimGalore
// =============================

process trim_galore {

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        path("*_val_1.fq.gz"), emit: trimmed_reads_R1
        path("*_val_2.fq.gz"), emit: trimmed_reads_R2

    publishDir params.trimgalore_output_dir, mode: params.publish_mode

    script:
    """
    echo "[DEBUG] Running Trim Galore on: $read1 and $read2 for sample $sample_id"
    trim_galore --paired "$read1" "$read2" --output_dir .
    """
}
// =============================
// PROCESS: PICARD FastqToSam (paired-end)
// =============================
process fastqtosam {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path("${sample_id}.bam")

    publishDir params.fastqtosam_output_dir, mode: params.publish_mode

    script:
    """
    picard FastqToSam \\
        F1=${read1} \\
        F2=${read2} \\
        O=${sample_id}.bam \\
        SM=${sample_id} \\
        SORT_ORDER=queryname
    """
}
