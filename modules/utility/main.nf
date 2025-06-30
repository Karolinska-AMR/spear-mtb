process GENERATE_REPORT{
  
    label 'process_low'

    input:
     path pip_out
     val prefix

    publishDir "${params.out_dir}/results", mode: 'copy', pattern: "${prefix}*.json"
    
    output:
     path("${prefix}*.json"), emit: json
    
    script:
    """
     python ${baseDir}/bin/generate_report.py --in_dir . --prefix $prefix
    """
}

process ARCHIVE_RAW_SEQ{
    tag "$sample"

    maxForks 1

    input:
    tuple val(sample), path(reads)
    val out_dir
    val cmd

    script:
    def dest= "${out_dir}/${sample}"
    """
    mkdir -p $dest
    r1=\$(readlink ${reads[0]})
    r2=\$(readlink ${reads[1]})
    ${cmd} \$r1 ${dest}/${sample}_1.fastq.gz
    ${cmd} \$r2 ${dest}/${sample}_2.fastq.gz
    """
}