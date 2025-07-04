process GENERATE_REPORT{
  
    input:
     path pip_out
     path html

    output:
    path("spear-mtb*.html"), emit: html
    
    script:
    """
     python ${baseDir}/bin/generate_report.py --in_dir . --html $html
    """
}