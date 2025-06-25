nextflow.enable.dsl = 2


// =============================
// WORKFLOW
// =============================
workflow {

    // ------------------------------------------
    // STEP 1: CONCATENATE RAW FASTQ FILES
    // ------------------------------------------
    raw_reads_ch = Channel.fromPath("${params.raw_fastq_input_dir}/*.fastq.gz")
        .filter { file ->
            def name = file.getName()
            return name ==~ /^HTL.*\.fastq\.gz$/
        }
        .map { file ->
            def name = file.getBaseName()
            // Removes the ".fastq.gz" suffix
            def tokens = name.split('_')
            // Split on underscores
            def sample = tokens[0]
            // First token = sample name
            def read = tokens.find { it ==~ /^R[12]$/ } ?: 'R?'
            // Read pair
            tuple("${sample}_${read}", file)
        }
        .groupTuple(by: 0)
        
        .map { sample_read, files ->
        def logFile = new File("Logs/concat_input_preview.log")
        if (!logFile.exists() || logFile.length() == 0) {
            // Write timestamp only once at the top if file is empty or doesn't exist
            logFile.text = "========== ${new Date().format('yyyy-MM-dd HH:mm:ss')} ==========\n\n"
        }
        logFile.withWriterAppend { writer ->
            writer << "${sample_read}:\n"
            writer << files*.name.join('\n') + "\n\n"
        }
        tuple(sample_read, files)
    }

    // Step 5: Print output (for debugging/inspection)
    // .view { group ->
    //     "Grouped: ${group[0]} -> ${group[1]*.getName()}"
    // }


    // Step 5: Concatenate files for each sample+read group
    concat_fastq_out_ch = concat_fastq(raw_reads_ch)
        .map { file ->  // `file` is a Path object returned by the process
        def logFile = new File("Logs/fastq_input_preview.log")
         if (!logFile.exists() || logFile.length() == 0) {
            // Write timestamp once at the top of the log file
            logFile.text = "========== ${new Date().format('yyyy-MM-dd HH:mm:ss')} ==========\n\n"
        }
        logFile.withWriterAppend { writer ->
            writer << "${file.name}\n\n"
        }
        file
    }



    
    // ------------------------------------------
    // STEP 2: FASTQC ON CONCATENATED FILES
    // ------------------------------------------

    raw_fastqc_out_ch = fastqc(concat_fastq_out_ch)


    // ------------------------------------------
    // STEP 3: MULTIQC ON RAW FASTQ FILES
    // ------------------------------------------

    raw_multiqc_input_ch = raw_fastqc_out_ch.fastqc_htmls // Collect all HTML reports
    .mix(raw_fastqc_out_ch.fastqc_zips)     // Mix with ZIP reports for MultiQC
    .collect() // Collect into a single channel for MultiQC

    multiqc_raw(raw_multiqc_input_ch)

    // ==============================================
    // STEP 4: QC ON CONCATENATED FASTQ FILES WITH TRIM GALORE
    // ==============================================

    trim_galore_input_ch = Channel
        .fromFilePairs("${params.concat_fastq_dir}/*_R{1,2}.fastq.gz", flat: false)
        .filter { sample_id, reads -> reads.size() == 2 }
        .map { sample_id, reads -> 
            tuple(sample_id, reads[0], reads[1])
        }

        // for debugging, print the input to Trim Galore
        // .view { "TrimGalore input tuple: ${it[0]}, R1: ${it[1].getName()}, R2: ${it[2].getName()}" }

    trim_galore_output_ch = trim_galore(trim_galore_input_ch)

    // ------------------------------------------
    // STEP 5: MULTIQC ON TRIMMED FILES
    // ------------------------------------------

   trimmed_fastqc_ch = trim_galore_output_ch.fastqc_htmls
    // Get the FastQC HTML outputs from Trim Galore

        .mix(trim_galore_output_ch.fastqc_zips)
        // Combine them with the FastQC ZIP outputs
        .mix(trim_galore_output_ch.trimming_reports)
        // Mix in the trimming reports for MultiQC
        .mix(trim_galore_output_ch.trimmed_reads_R1)
        // Include the trimmed reads R1         
        .mix(trim_galore_output_ch.trimmed_reads_R2)
        // Include the trimmed reads R2
        .flatten()
        // Convert grouped file lists into a flat stream of individual files
    
        .collect()
        // Gather all individual files into a single list for MultiQC
    
        // .view { "Trimmed FastQC output:\n" + it*.getName().join('\n') }
        // (Optional) Print all collected filenames for debugging

    multiqc_trimmed(trimmed_fastqc_ch)
    // Run MultiQC on the collected FastQC outputs


    // ------------------------------------------
    // STEP 4: PICARD FastqToSam (paired-end)
    // ------------------------------------------
    // fastqtosam_input_ch = Channel
    // .fromFilePairs("${params.trimgalore_output_dir}/HTL284*R{1,2}_val_{1,2}.fq.gz", flat: true)
    // .map { sample_id, reads -> 
    //     tuple(sample_id, reads[0], reads[1])
    // }
    // .view()
    
  

    // fastqtosam(fastqtosam_input_ch)
}








// =============================
// PROCESS: CONCATENATE FASTQ
// =============================
process concat_fastq {

    tag "${sample_read}"

    publishDir params.concat_fastq_dir, mode: params.publish_mode

    input:
    tuple val(sample_read), path(fastq_files)

    output:
    path "${sample_read}.fastq.gz"

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

    publishDir params.fastqc_dir, mode: params.publish_mode

    input:
    path input_file

    output:
    path "*fastqc.html", emit: fastqc_htmls
    path "*fastqc.zip", emit: fastqc_zips

    script:
    """
    echo "🔍 Running FastQC on ${input_file.getName()}" >&2
    fastqc ${input_file} --outdir .
    echo "✅ FastQC done for ${input_file.getName()}" >&2
    """
}

// ====================================
// PROCESS: MULTIQC ON RAW FASTQC FILES
// ====================================

process multiqc_raw {

    tag "multiqc_raw_fastq"

    publishDir params.multiqc_raw_dir, mode: params.publish_mode

    input:
    path fastqc_dirs

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . -o . 
    """
}

// ========================================
// PROCESS: MULTIQC ON TRIMMED FASTQC FILES
// ========================================

process multiqc_trimmed {

    tag "multiqc_trimmed_fastq"

    publishDir params.multiqc_trimmed_dir, mode: params.publish_mode

    input:
    path fastqc_dirs

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . -o . 
    """
}

// =============================
// PROCESS: TrimGalore
// =============================

process trim_galore {

    tag "${sample_id}"

    publishDir "QC/trimmed_fastqc_reports", pattern: "*_fastqc.*", mode: params.publish_mode
    publishDir "QC/trimmed_fastqc_reports", pattern: "*trimming_report.txt", mode: params.publish_mode
    publishDir params.trimgalore_output_dir, pattern: "*_val_*.fq.gz", mode: params.publish_mode

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path ("*_val_1.fq.gz"), emit: trimmed_reads_R1
    path ("*_val_2.fq.gz"), emit: trimmed_reads_R2
    path ("*_fastqc.html"), emit: fastqc_htmls
    path ("*_fastqc.zip"), emit: fastqc_zips
    path ("*trimming_report.txt"), optional: true, emit: trimming_reports
    path ("trim_galore_${sample_id}.log"), emit: stdout_log
    path ("trim_galore_${sample_id}.err"), emit: stderr_log

    script:
    """
    echo "[DEBUG] Running Trim Galore on: ${read1} and ${read2} for sample ${sample_id}"

    trim_galore \\
        --paired \\
        -a " AGATCGGAAGAGC -a G{50}" \\
        -a2 " AGATCGGAAGAGC -a G{50}" \\
        --quality 20 \\
        --clip_R1 5 \\
        --clip_R2 5 \\
        --length 20 \\
        --cores 4 \\
        --gzip \\
        --fastqc \\
        --output_dir . \\
        "${read1}" "${read2}" \\
        > trim_galore_${sample_id}.log 2> trim_galore_${sample_id}.err
    """
}
// =============================
// PROCESS: PICARD FastqToSam (paired-end)
// =============================
process fastqtosam {

    tag { "${sample_id}" }

    publishDir params.fastqtosam_output_dir, mode: params.publish_mode

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "${sample_id}.bam"

    script:
    """
    echo "[DEBUG] Running Picard FastqToSam for sample: ${sample_id}"

    picard FastqToSam \\
        F1=${read1} \\
        F2=${read2} \\
        O=${sample_id}_unmapped.bam \\
        SM=${sample_id} \\
        RG=${sample_id} \\
        PL=ILLUMINA \\
        SORT_ORDER=queryname \\
        REFERENCE_SEQUENCE="${params.reference_genome}" \\
        TMP_DIR="${params.tmp_dir}"
    """

}
