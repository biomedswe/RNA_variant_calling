nextflow.enable.dsl=2

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
            // Accept files that match: HTL123_Sx_Lxxx_R[12]_xxx.fastq.gz
            return name ==~ /^HTL\d+.*_R[12]_.*\.fastq\.gz$/
        }

        // Step 2: Extract sample and read (R1 or R2)
        .map { file ->
            def name = file.getName()
            def tokens = name.split('_')
            def sample = tokens[0]
            def read   = tokens[3]
            tuple("${sample}_${read}", file)
        }

        // Step 3: Group files by sample+read (e.g. HTL123_R1)
        .groupTuple()

        // Step 4: Concatenate grouped FASTQ files
        | concat_fastq
}

process concat_fastq {

    // Limit how many concat_fastq jobs run at the same time (avoid overloading cluster)
    maxForks 20

    // Resource requests per job
    cpus 1                      // Number of CPU cores per job
    memory '2 GB'               // Amount of memory requested
    time '1h'                   // Maximum walltime per job (if scheduler enforces it)

    // 🏷 Tag for easier tracking in logs
    tag "${group}"

    // Inputs: a tuple of sample_read group and a list of FASTQ files
    input:
    tuple val(group), path(fastq_files)

    // Output: the final merged FASTQ file
    output:
    path("${group}.fastq.gz")

    // Copy output to target folder (rather than leaving it in work/)
    publishDir params.output_dir, mode: 'copy'

    

    // Skip this job if its output file already exists
    when:
    !file("${params.output_dir}/${group}.fastq.gz").exists()

    // The actual shell script to run
    script:
    """
    mkdir -p ${params.log_dir}

    {
    cat ${fastq_files.join(' ')} > ${group}.fastq.gz
    } > ${params.log_dir}/${group}.out 2> ${params.log_dir}/${group}.err
    """
}
