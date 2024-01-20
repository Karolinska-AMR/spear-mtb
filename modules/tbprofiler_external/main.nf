process TBPROFILER_PROFILE_EXTERNAL {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tb-profiler:5.0.1--pyhdfd78af_1' :
        'biocontainers/tb-profiler:5.0.1--pyhdfd78af_1' }"
    containerOptions "--bind $db:/usr/local/share/tbprofiler"

    input:
    tuple val(meta), path(bam), path(db)

    output:
    tuple val(meta), path("results/*.csv") , emit: csv, optional: true
    tuple val(meta), path("results/*.json"), emit: json
    tuple val(meta), path("results/*.txt") , emit: txt, optional: true
    path "versions.yml"                    , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tbprofiler:  \$( echo \$(tb-profiler --version 2>&1) | sed 's/TBProfiler version //')
    END_VERSIONS
    """
}