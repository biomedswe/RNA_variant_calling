// =============================
// PROCESS: FASTQC
// =============================
process FastQC {

    tag { sample_read }

    publishDir { "QC/${output_subfolder}/${sample}" }, mode: params.publish_mode

    input:
    tuple val(sample), val(sample_read), path(read)
    val(output_subfolder)

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
    
    if ! fastqc ${read} --outdir . 1>> \$PWD/.command.out 2>&1; then
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