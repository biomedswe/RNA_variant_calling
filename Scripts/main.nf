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



    
    // // ------------------------------------------
    // // STEP 2: FASTQC ON CONCATENATED FILES
    // // ------------------------------------------

    raw_fastqc_out_ch = fastqc(concat_fastq_out_ch)


    // // ------------------------------------------
    // // STEP 3: MULTIQC ON RAW FASTQ FILES
    // // ------------------------------------------

    // // print the collected FastQC outputs for debugging
    // // .collect consume the channel and return a single list of all FastQC outputs
    // multiqc_input_ch = raw_fastqc_out_ch.collect()
    
    // // multiqc_input_ch.view { collected ->
    // // "Collected FastQC outputs:\n" + collected*.getName().join("\n")
    // // }

    // // Collect all outputs from fastqc
    // multiqc(multiqc_input_ch) // Collects all FastQC outputs and sends them to MultiQC. without .collect() it would send each file separately, which is not what we want.

    // // multiqc_input_ch = fastqc_out_ch.flatten().collect()
    // // multiqc_input_ch.view { collected ->
    // // "Collected FastQC outputs:\n" + collected*.getName().join("\n")
    // // }
    // // multiqc(multiqc_input_ch)

    raw_multiqc_input_ch = raw_fastqc_out_ch.fastqc_htmls
    .mix(raw_fastqc_out_ch.fastqc_zips)
    .collect()

    multiqc(raw_multiqc_input_ch)

    // // ==============================================
    // // STEP 4: QC ON CONCATENATED FASTQ FILES WITH TRIM GALORE
    // // ==============================================

    // trim_galore_input_ch = Channel
    //     .fromFilePairs("${params.concat_fastq_dir}/*_R{1,2}.fastq.gz", flat: false)
    //     .filter { sample_id, reads -> reads.size() == 2 }
    //     .map { sample_id, reads -> 
    //         tuple(sample_id, reads[0], reads[1])
    //     }

    //     // for debugging, print the input to Trim Galore
    //     // .view { "TrimGalore input tuple: ${it[0]}, R1: ${it[1].getName()}, R2: ${it[2].getName()}" }

    // trim_galore_output_ch = trim_galore(trim_galore_input_ch)

    
    // ==============================================
    // STEP 4: QC ON RAW FASTQ FILES WITH TRIM GALORE
    // ==============================================

// trim_galore_input_ch = concat_fastq_out_ch
//     // Step 1: Strip _R1/_R2 to get sample ID
//     .map { sample_read, file -> 
//         def sample_id = sample_read.replaceAll(/_R[12]$/, '')
//         tuple(sample_id, file)
//     }
//     .view { "After .map → (sample_id, file): $it" }

    // // Step 2: Group all files for the same sample
    // .groupTuple(by: 0)
    // .view { "After .groupTuple → (sample_id, [files]): $it" }

    // // Step 3: Find R1 and R2 files and prepare for TrimGalore
    // .map { sample_id, files ->
    //     def r1 = files.find { it.getFileName().toString().contains('_R1') }
    //     def r2 = files.find { it.getFileName().toString().contains('_R2') }
    //     tuple(sample_id, r1, r2)
    // }
    // .view { "Final TrimGalore input → (sample_id, R1, R2): $it" }



    // Start with the output from concat_fastq (one merged FASTQ file per R1 or R2 read)
    // Example input: ("HTL123_R1", path/to/HTL123_R1.fastq.gz)
    // trim_galore_input_ch = concat_fastq_out_ch
    // // Converts: ("HTL123_R1", file1) → ("HTL123", file1)
    // // Converts: ("HTL123_R2", file2) → ("HTL123", file2)
    // // This strips the _R1 or _R2 from the filename to extract a common sample ID (like "HTL123"), so both files can be grouped together.
    // .map { sample_read, file -> 
    //     def sample_id = sample_read.replaceAll(/_R[12]$/, '')  // Strip _R1 or _R2
    //     tuple(sample_id, file)  // Return a tuple: (sample_id, file)
    // }

    // .groupTuple(by: 0)
   
    // // Example: ("HTL123", [("HTL123", file1), ("HTL123", file2)]) → ("HTL123", file1, file2)
    // .map { sample_read, file_list ->
    // def r1 = file_list.find { it[1].name.contains('_R1') }[1]
    // def r2 = file_list.find { it[1].name.contains('_R2') }[1]
    // tuple(sample_read, r1, r2)
    // }

    // // // This final channel (trim_input_ch) looks like:
    // // // ("HTL123", file1_R1.fastq.gz, file2_R2.fastq.gz)
    // // // and is sent into the trim_galore process
    // trim_galore_output_ch = trim_galore(trim_galore_input_ch)

    // // ------------------------------------------
    // // STEP 4: MULTIQC ON TRIMMED FILES
    // // ------------------------------------------

    // // mix combines the two channels into one, so MultiQC can process both HTML and ZIP files together
    // trimmed_fastqc_ch = trim_galore_output_ch.fastqc_htmls
    //     .mix(trim_galore_output_ch.fastqc_zips)
    //     .view { "Trimmed FastQC output: ${it.getName()}" }

    // multiqc(trimmed_fastqc_ch)
  

    

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

// =============================
// PROCESS: MULTIQC
// =============================

process multiqc {

    tag "multiqc_raw_fastq"

    publishDir params.multiqc_dir, mode: params.publish_mode

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

    # Normalize input filenames using symlinks
    ln -s "${read1}" "${sample_id}_1.fastq.gz"
    ln -s "${read2}" "${sample_id}_2.fastq.gz"

    # Run Trim Galore with FastQC
    # Capture stdout and stderr
    trim_galore --paired "${sample_id}_1.fastq.gz" "${sample_id}_2.fastq.gz" --gzip --fastqc --output_dir . \
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
    picard FastqToSam \\
        F1=${read1} \\
        F2=${read2} \\
        O=${sample_id}.bam \\
        SM=${sample_id} \\
        SORT_ORDER=queryname
    """
}
