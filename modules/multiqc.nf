/*
 * modules/multiqc.nf
 * A single-purpose, re-usable MultiQC process
 */
 
process MultiQC {

    tag { "${report_tag}" }

    publishDir "QC/${report_tag}", mode: params.publish_mode

    input:
    tuple val(report_tag), path(fastqc_dirs)

    output:
    path "multiqc_report.html"

    

    script:
    """
    # Source shared functions
    source ${params.script_dir}/nextflow_functions.sh

    # Use log_dir_name from workflow
    export LOG_DIR="${params.log_dir}/${report_tag}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "🔍 Running ${task.process} on ${task.tag}" >> "\$PWD/.command.out"

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