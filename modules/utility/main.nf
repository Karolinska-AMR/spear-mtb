process GENERATE_REPORT{
  
    input:
     path pip_out
     val prefix
    output:
     path("${prefix}*.json"), emit: json
    script:
    """
     python ${baseDir}/bin/generate_report.py --in_dir .  --prefix $prefix
    """
}