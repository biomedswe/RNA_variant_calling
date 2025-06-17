nextflow.enable.dsl = 2

// Create missing output/log folders BEFORE anything runs
// Runs in Nextflow's Groovy context (not as a job)
[
    params.output_dir,
    params.log_dir
].each { dir ->
    def path = file(dir)
    if (!path.exists()) {
        println "Creating missing directory: ${dir}"
        path.mkdirs()
    }
}

workflow {
    Channel
        .fromPath("${params.input_dir}/*.fastq.gz")

        // Step 1: Validate filenames with expected pattern
        .filter { file ->
            def name = file.getName()
            return name ==~ /^HTL.*\.fastq\.gz$/ // Ensure the file name starts with "HTL" and ends with ".fastq.gz"
        }

// Step 1: Read FASTQ files and extract a clean sample + read name
    .map { file ->
            def name = file.getBaseName()  // Removes the ".fastq.gz" suffix
                                           // Example: "HTL284_S1_L001_R1_001.fastq.gz" → "HTL284_S1_L001_R1_001"

            def tokens = name.split('_')  // Split on underscores to break up the name into parts

            def sample = tokens[0]        // First token is typically the sample name, e.g., "HTL284"

            // Try to find "R1" or "R2" in the tokens; if not found, set to "R?" (helps spot broken names)
            def read = tokens.find { it ==~ /^R[12]$/ } ?: 'R?'

            // Combine sample and read into a grouping key, like "HTL284_R1"
            tuple("${sample}_${read}", file)
}

        // Step 3: Group files by sample+read (e.g. HTL123_R1)
        .groupTuple()

        // Step 4: Concatenate grouped FASTQ files
        .filter { group ->
            def sample_read = group[0]
            def output_path = file("${params.output_dir}/${sample_read}.fastq.gz")
            if (output_path.exists()) {
                println "⚠️  Skipping ${sample_read} — file already exists."
                return false
            }
            return true
        }

        | concat_fastq
}

process concat_fastq {

    // =============================
    // RESOURCE ALLOCATION
    // =============================
    maxForks 100        // Allow up to 100 concurrent jobs (tune this based on your machine/cluster)
    cpus 1              // Use 1 CPU per job
    memory '2 GB'       // Allocate 2 GB RAM
    time '1h'           // Set max runtime per job (useful on clusters)

    // =============================
    // PROCESS IDENTIFIER
    // =============================
    tag "${sample_read}"  // This will show in logs, trace reports, etc. Helpful for debugging

    // =============================
    // INPUTS: A tuple of (sample_read name, list of FASTQ files)
    // =============================
    input:
    tuple val(sample_read), path(fastq_files)

    // =============================
    // OUTPUT: The resulting merged file named like HTL123_R1.fastq.gz
    // =============================
    output:
    path("${sample_read}.fastq.gz")

    // =============================
    // TARGET OUTPUT FOLDER
    // =============================
    publishDir params.output_dir, mode: 'copy'

   

    // =============================
    // PROCESS BODY: merge or copy input FASTQ files
    // =============================
    script:
    def output_file = "${sample_read}.fastq.gz"

    // Build the shell command:
    // If only one file, just copy it
    // If multiple files (e.g. different lanes), concatenate them
    def cmd = fastq_files.size() == 1 
        ? "cp ${fastq_files[0]} ${output_file}" 
        : "cat ${fastq_files.join(' ')} > ${output_file}"

    // Run the command
    // Use >&2 to send log/debug info to stderr (Nextflow shows it in terminal)
    """
    mkdir -p ${params.log_dir}

    echo "→ Merging files for ${sample_read}" >&2
    ${cmd}
    echo "✅ Created: ${output_file}" >&2
    """
}