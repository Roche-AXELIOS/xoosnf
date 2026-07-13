process XOOS_DEMUX {
    conda false
    container 'ghcr.io/roche-axelios/xoos/demux:1.1.0'

    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(fastq_dir), path(sample_sheet)
    path strand_detect_index, stageAs: 'strand_index/*'

    output:
    tuple val(meta), path('demux/*'), emit: demux_output
    tuple val(meta), path('demux/metrics'), emit: demux_metrics
    tuple val("${task.process}"), val('demux'), eval('demux --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def strand_detect_args = strand_detect_index ? "--strand-detect ${strand_detect_index}" : ''
    def args = task.ext.args ?: ''
    """
    demux \\
        --input '${fastq_dir}' \\
        --out-dir 'demux' \\
        --sample-sheet '${sample_sheet}' \\
        --adapter-design-name '${meta.run_type}' \\
        ${strand_detect_args} \\
        ${args}
    """

    stub:
    def stub_args = task.ext.args ?: ''
    def is_all_split = stub_args.contains('all-split')
    """
    install -D /dev/null demux/metrics/run_stats.tsv
    install -D /dev/null demux/metrics/sample_stats.tsv
    tail -n +2 ${sample_sheet} | while IFS=, read -r sample_name sample_sid || [[ -n "\$sample_name" ]]; do
        # is_all_split is interpolated as 'true' or 'false' (valid bash builtins)
        if ${is_all_split}; then
            mkdir -p "demux/\${sample_name}/full"
            mkdir -p "demux/\${sample_name}/partial"
            touch demux/\${sample_name}/full/\${sample_name}_001.fastq.gz
            touch demux/\${sample_name}/full/\${sample_name}_002.fastq.gz
            touch demux/\${sample_name}/partial/\${sample_name}_001.fastq.gz
        else
            mkdir -p "demux/\${sample_name}"
            touch demux/\${sample_name}/\${sample_name}_001.fastq.gz
            touch demux/\${sample_name}/\${sample_name}_002.fastq.gz
        fi
    done
    """
}
