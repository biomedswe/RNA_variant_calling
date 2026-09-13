// =============================
// PROCESS: FeatureCounts
// =============================

process FeatureCounts {

    tag { "${strandedness}" }

    input:
    path(wasp_sort_dedup_folder)
    val(strandedness)

    output:
    tuple val(strandedness), path("HTL2_FeatureCounts_${strandedness}.txt")

    script:
    """
    # Source shared functions
    source ${file("${moduleDir}/nextflow_functions.sh")}

    export LOG_DIR="${params.log_dir}/${task.process}"
    export TAG="${task.tag}"
    mkdir -p "\$LOG_DIR"

    echo "[DEBUG] 🔁 ${task.process} bam for sample: ${task.tag}" >> "\$PWD/.command.out"

    # run FeatureCounts
    if ! featureCounts \
        -a ${params.reference_gtf} \
        -o HTL2_FeatureCounts_${strandedness}.txt \
        ${wasp_sort_dedup_folder} \
        -p -T 12 -s ${strandedness} -C -O \
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