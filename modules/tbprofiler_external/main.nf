process TBPROFILER_PROFILE_EXTERNAL {
    tag "$meta.id"
    label 'process_medium'
    memory { 5.GB * task.attempt }
	maxRetries 3
	errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:5.0.1--pyhdfd78af_1' :
        'biocontainers/tb-profiler:5.0.1--pyhdfd78af_1' }"
    containerOptions "--bind $db:/usr/local/share/tbprofiler"

    input:
    tuple val(meta), path(bam), path(db)
    val pub_dir

    output:
    tuple val(meta), path("results/*.json"), emit: json
  

    publishDir { "${pub_dir}/${meta.id.tokenize('.')[0]}/tbprofiler"},
            mode: 'copy',
            saveAs: { fn -> fn.tokenize('/')[-1] }
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"    
    
    """   
    
    db_pfx=`find -L ./$db -name "*.fasta" | xargs -I{} basename {} .fasta`
 
    tb-profiler \\
        profile \\
        $args \\
        --prefix ${prefix} \\
        --threads $task.cpus \\
        --bam $bam \\
        --db \$db_pfx

    """
}