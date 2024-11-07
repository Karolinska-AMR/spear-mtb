nextflow.enable.dsl=2
include {KRAKEN2_KRAKEN2} from "./modules/kraken2"
include {BRACKEN_BRACKEN} from "./modules/braken"



process MERGE_BARACKEN{
    label 'process_low'

    input:
    path tsvs
    

    output:
    path "species_samples.out", emit: tsv
    
    
    script:
    
    
    
    """
    result_file="species_samples.out"
    touch "\$result_file"
    header_added=false

    # Iterate through each file in the directory
    for file in ./*.tsv; do
        # Get the file name without the path
        filename=\$(basename "\$file")
        sample="\${filename%.b*}"

        if [ "\$header_added" = false ]; then
                header=\$(head -n1 \$file)
                echo -e "\$header\tsample" > "\$result_file"
                header_added=true
        fi
        # Sort the file based on the last column and retrieve the first row
        first_row=\$(sort -t\$'\t' -k7,7r "\$file" | head -n 2 | tail -n 1 )
        
        # Append the file name as another column
        echo -e "\$first_row\t\$sample" >> "\$result_file"

    done
    """
    
}

workflow MTB_FINDER{
    
    take:
        input_ch
        kraken2_db
   
    main:
        input_ch = input_ch.map{it->[[id:it[0]],it[1]]}                     
    //Running Kraken2
        KRAKEN2_KRAKEN2(input_ch,kraken2_db,false,false)

        // Defining taxonomy levels
        level_ch = Channel.of('S')
        bracken_ch = KRAKEN2_KRAKEN2.out.report.combine(level_ch)
        
        // Running braken
        BRACKEN_BRACKEN(bracken_ch,kraken2_db)
        tsv_ch = BRACKEN_BRACKEN.out.reports.map{it->it[1]}.flatten().collect()   

        MERGE_BARACKEN(tsv_ch)
        mtbs_smaples = MERGE_BARACKEN.out.tsv.splitCsv(sep:"\t",header:true)
                    .filter{row -> row.fraction_total_reads.toDouble() >= 0.8 && row.name == "Mycobacterium tuberculosis"}
                    .map{it -> [[id:it.sample]]}
        
        mtbs_ch = input_ch.join(mtbs_smaples).map{it->[it[0].id,it[1]]}

    emit:
        species_tsvs = MERGE_BARACKEN.out.tsv
        mtbs_ch = mtbs_ch
}
