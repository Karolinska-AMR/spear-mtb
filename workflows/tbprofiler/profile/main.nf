process TBPROFILER_PROFILE {
    tag "$meta.id"
    label 'process_medium'
	 
	memory { 5.GB * task.attempt }
	maxRetries 2
	errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
	 
    conda (params.enable_conda ? "bioconda::tb-profiler=5.0.1-pyhdfd78af_1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:5.0.1--pyhdfd78af_1' :
        'quay.io/biocontainers/tb-profiler:5.0.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(reads)
    val pub_dir

    output:
    tuple val(meta), path("bam/*.bam")     , emit: bam
    tuple val(meta), path("bam/*.bam.bai") , emit: bai
    tuple val(meta), path("results/*.csv") , emit: csv, optional: true
    tuple val(meta), path("results/*.json"), emit: json
    tuple val(meta), path("results/*.txt") , emit: txt, optional: true
    tuple val(meta), path("vcf/*.vcf.gz")  , emit: vcf
    path "versions.yml"                    , emit: versions

    publishDir { "${pub_dir}/${meta.id.tokenize('.')[0]}/tbprofiler"},
            mode: 'copy',
            saveAs: { fn -> fn.tokenize('/')[-1] }


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "--read1 $reads" : "--read1 ${reads[0]} --read2 ${reads[1]}"
    """
    tb-profiler \\
        profile \\
        $args \\
        --prefix ${prefix} \\
        --threads $task.cpus \\
        $input_reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tbprofiler:  \$( echo \$(tb-profiler --version 2>&1) | sed 's/TBProfiler version //')
    END_VERSIONS
    """
}