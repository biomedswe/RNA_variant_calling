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
        .filter { sample_read, files -> sample_read.startsWith('HTL214') } // <-- Filter for HTL214
        
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
    concat_fastq_out_ch = Concatenating_fastq(raw_reads_ch)
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

    raw_fastqc_out_ch = FastQC(concat_fastq_out_ch)


    // ------------------------------------------
    // STEP 3: MULTIQC ON RAW FASTQ FILES
    // ------------------------------------------

    raw_multiqc_input_ch = raw_fastqc_out_ch.fastqc_htmls // Collect all HTML reports
    .mix(raw_fastqc_out_ch.fastqc_zips)     // Mix with ZIP reports for MultiQC
    .collect() // Collect into a single channel for MultiQC

    multiqc_raw(raw_multiqc_input_ch)

    // =======================================================
    // STEP 4: QC ON CONCATENATED FASTQ FILES WITH TRIM GALORE
    // =======================================================

    // trim_galore_input_ch = Channel
    //     .fromFilePairs("${params.concat_fastq_dir}/*_R{1,2}.fastq.gz", flat: false)
    //     .filter { sample_id, reads -> reads.size() == 2 }
    //     .map { sample_id, reads -> 
    //         tuple(sample_id, reads[0], reads[1])
    //     }

    trim_galore_input_ch = Channel
    // 1. collect only proper R1/R2 pairs
    .fromFilePairs("${params.concat_fastq_dir}/*_R{1,2}.fastq.gz")
    // 2. keep only the samples of interest
    .filter { sample_id, reads -> sample_id.startsWith('HTL214') }
    // 3. explode the pair into three elements if that’s what the next
    //    process expects (otherwise drop this map)
    .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
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
    // STEP 6: STAR ALIGNMENT
    // ------------------------------------------

    // Join the R1 and R2 files on their sample prefix
    fastqtosam_align_reads_input_ch = trim_galore_output_ch.trimmed_reads_R1
        .combine(trim_galore_output_ch.trimmed_reads_R2)
        .map { r1, r2 ->
            def sample_id = r1.getName().replaceFirst(/_R1_val_1\.fq\.gz$/, '') // Full file name with extensions
            tuple(sample_id, r1, r2)
        }
        
    // Input trimmed fastq files from trim_galore
    alignreads_output_ch = AlignReads(fastqtosam_align_reads_input_ch)

    // View AlignReads output
    // alignreads_ch.view { ">> alignreads_ch: ${it}" }

    // ------------------------------------------
    // STEP 7: MarkDuplicates
    // ------------------------------------------

    splitncigarreads_input_ch = MarkDuplicates(alignreads_output_ch)
    .map { sample_id, bam_file, bai_file ->
        tuple(sample_id, bam_file)  // drop the bai_file
    }

    
    // ------------------------------------------
    // STEP 8: SplitNCigarReads
    // ------------------------------------------

    splitncigarreads_output_ch = SplitNCigarReads(splitncigarreads_input_ch)
    splitncigarreads_output_ch.view { "splitncigarreads_output_ch: ${it}"}

    // ------------------------------------------
    // STEP 9: BaseRecalibrator
    // ------------------------------------------

    base_recalibrator_out_ch = BaseRecalibrator(splitncigarreads_output_ch)

    base_recalibrator_out_ch.view { "base_recalibrator_out_ch: ${it}"}

    // Remap to pass BAM + recal table to ApplyBQSR
//    apply_bqsr_input_ch = base_recalibrator_out_ch
//     .join(splitncigarreads_output_ch) // .join is best practice for joining channels when they are related
//     .map { sample_id, recal, bam ->
//         def recal_table = recal[1]
//         def bam_file    = bam[1]
//         tuple(sample_id, bam_file, recal_table)
//     }
//     .view { "apply_bqsr_input_ch: ${it}" }

    apply_bqsr_input_ch = base_recalibrator_out_ch
        .join(splitncigarreads_output_ch)
        .map { sample_id, recal_table, bam_file ->
            tuple(sample_id, bam_file, recal_table)
        }
        .view { "apply_bqsr_input_ch: ${it}" }

// Step 10: ApplyBQSR
ApplyBQSR(apply_bqsr_input_ch)

}









// =============================
// PROCESS: CONCATENATE FASTQ
// =============================
process Concatenating_fastq {

    tag "${sample_read}"

    input:
    tuple val(sample_read), path(fastq_files)

    output:
    path "${sample_read}.fastq.gz"

    script:
    // This is groovy code and therefore outside the triple quotes
    // Define the output file name based on the sample name
    def output_file = "${sample_read}.fastq.gz"
       
    // - If there's only one input FASTQ file, just copy it (no need to concatenate).
    // - If there are multiple FASTQ files (e.g., multiple lanes for the same sample),
    //   concatenate them into a single output file using 'cat'.
    def cmd = fastq_files.size() == 1
        ? "cp ${fastq_files[0]} ${output_file}"                        // Case: only one file — use 'cp'
        : "cat ${fastq_files.join(' ')} > ${output_file}"              // Case: multiple files — use 'cat' to merge

    

    """
    # Source shared functions
    source ${params.script_dir}/nextflow_functions.sh

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
// PROCESS: FASTQC
// =============================
process FastQC {

    tag { input_file.getBaseName() }

    input:
    path input_file

    output:
    path "*fastqc.html", emit: fastqc_htmls
    path "*fastqc.zip", emit: fastqc_zips

    script:

    
    """
    # Source shared functions
    source ${params.script_dir}/nextflow_functions.sh

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "🔍 Running ${task.process} on ${task.tag}" >> "\$PWD/.command.out"
    
    if ! fastqc ${input_file} --outdir . >> \$PWD/.command.out 2>> \$PWD/.command.err ; then
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} was succesful for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
    """
}

// ====================================
// PROCESS: MULTIQC ON RAW FASTQC FILES
// ====================================

process multiqc_raw {

    tag "multiqc_raw_fastq"

    input:
    path fastqc_dirs

    output:
    path "multiqc_report.html"

    script:
    

    """
    # Source shared functions
    source ${params.script_dir}/nextflow_functions.sh

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "🔍 Running ${task.process} on ${task.tag}" >> "\$PWD/.command.out"
    
    # append stdout and stderr
    if ! multiqc . -o . 1>> \$PWD/.command.out 2>&1; then
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} was succesful for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
     
    """
}

// ========================================
// PROCESS: MULTIQC ON TRIMMED FASTQC FILES
// ========================================

process multiqc_trimmed {

    tag "multiqc_trimmed_fastq"

    input:
    path fastqc_dirs

    output:
    path "multiqc_report.html"

    script:
    

    """
    # Source shared functions
    source ${params.script_dir}/nextflow_functions.sh

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "🔍 Running ${task.process} on ${task.tag}" >> "\$PWD/.command.out"
    
    # append stdout and stderr
    if ! multiqc . -o . 1>> \$PWD/.command.out 2>&1; then
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    else
        echo "✅ ${task.process} was succesful for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0
    fi
     
    """
}

// =============================
// PROCESS: TrimGalore
// =============================

process trim_galore {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path ("*_val_1.fq.gz"), emit: trimmed_reads_R1
    path ("*_val_2.fq.gz"), emit: trimmed_reads_R2
    path ("*_fastqc.html"), emit: fastqc_htmls
    path ("*_fastqc.zip"), emit: fastqc_zips
    path ("*trimming_report.txt"), optional: true, emit: trimming_reports
    

    script:
     
    """
    # Source shared functions
    source ${params.script_dir}/nextflow_functions.sh

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"
   
    echo "🔍 Running ${task.process} on: ${read1} and ${read2} for sample ${task.tag}"

    if ! trim_galore \
        --paired \
        -a " AGATCGGAAGAGC -a G{50}" \
        -a2 " AGATCGGAAGAGC -a G{50}" \
        --quality 20 \
        --clip_R1 5 \
        --clip_R2 5 \
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


// ==========================
// PROCESS: STAR alignReads 
// ==========================


process AlignReads {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    

    """
    # Source shared functions
    source ${params.script_dir}/nextflow_functions.sh

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
        outFilterMismatchNoverLmax 0.04 \
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
    source ${params.script_dir}/nextflow_functions.sh

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
    source ${params.script_dir}/nextflow_functions.sh

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
    source ${params.script_dir}/nextflow_functions.sh

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
    source ${params.script_dir}/nextflow_functions.sh

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
    source ${params.script_dir}/nextflow_functions.sh

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
        -jar /ludc/Home/jonas_a/.conda/envs/variant_calling/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
        HaplotypeCaller \
        -R ${params.reference_genome} \
        -I ${input_bam} \
        -O "${sample_id}_haplotypecaller.vcf.gz" \
        --dont-use-soft-clipped-bases \
        --dbsnp ${params.snp_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
        1>> \$PWD/.command.out 2>&1; then

        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0

    else
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
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
    source ${params.script_dir}/nextflow_functions.sh

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Run GATK VariantFiltration
    if ! java -Djava.io.tmpdir=${params.tmp_dir} \
        -Xms6000m \
        -XX:GCTimeLimit=50 \
        -XX:GCHeapFreeLimit=10 \
        -XX:+PrintFlagsFinal \
        -Xlog:gc*:file=gc_log.log:time,uptime,level,tags \
        -jar /ludc/Home/jonas_a/.conda/envs/variant_calling/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar \
        VariantFiltration \
        -R ${params.reference_genome} \
        -V ${input_vcf} \
        -window 35 \
        -cluster 3 \
        --filter-name "FS" -filter "FS > 30.0" \
        --filter-name "QD" -filter "QD < 2.0" \
        -O "${sample_id}_filtered.vcf.gz" \
        1>> \$PWD/.command.out 2>&1; then

        echo "✅ ${task.process} done for ${task.tag}" >> "\$PWD/.command.out"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 0

    else
        echo "❌ ${task.process} failed for ${task.tag}" >> "\$PWD/.command.err"
        move_named_log ".command.condor" ".command.err" ".command.out" ".command.run" ".command.sh"
        exit 1
    fi
    """
}

// =============================
// PROCESS: FindIntersectingSNPs
// =============================

process FindIntersectingSNPs {

    tag { "${sample_id}" }

    input:
    tuple val(sample_id), path(input_bam), path(filtered_vcf)

    output:
    tuple val(sample_id), path("*_remap.bam"), path("*_remap.fq.gz"), path("*_to_remap.bam")

    script:
    """
    # Source shared functions
    source ${params.script_dir}/nextflow_functions.sh

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔍 Running ${task.process} for sample: ${task.tag}" >> "\$PWD/.command.out"

    # Run find_intersecting_snps.py
    if ! python ${params.wasp_dir}/mapping/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir . \
        --snp_dir ${filtered_vcf} \
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







	