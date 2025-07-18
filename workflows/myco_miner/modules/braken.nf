process BRACKEN_BRACKEN {
    tag "$meta.id"
    label 'process_low'


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bracken:2.7--py39hc16433a_0':
        'quay.io/biocontainers/bracken:2.7--py39hc16433a_0' }"

    input:
    tuple val(meta), path(kraken_report), val(level)
    path database

    output:
    tuple val(meta), path(bracken_report), emit: reports
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    bracken_report = "${prefix}.braken-${level}.tsv"
    def VERSION = '2.7'
    """
    bracken \\
        ${args} \\
        -d '${database}' \\
        -i '${kraken_report}' \\
        -o '${bracken_report}' \\
        -l '${level}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: ${VERSION}
    END_VERSIONS
    """
}