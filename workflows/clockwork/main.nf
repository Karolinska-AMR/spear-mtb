 process MAP_READS {
    tag "$meta.id"
    label 'process_medium'
    
    memory { 5.GB * task.attempt }
	maxRetries 3
	errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/iqbal-lab-org/clockwork/releases/download/v0.11.3/clockwork_v0.11.3.img':
        'ghcr.io/iqbal-lab-org/clockwork:latest' }"

    input:
    tuple val(meta), path(reads), val(ref_fa)
       

    output:
    tuple val(meta), path("*.sam"), emit: sam

    script:
     def threads = task.cpus
    """
     clockwork map_reads --threads $threads --unsorted_sam ${meta.id} $ref_fa ${meta.id}.sam ${reads[0]} ${reads[1]}
    """
}
process REMOVE_CONTAM {
    tag "$meta.id"
    label 'process_medium'
   
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/iqbal-lab-org/clockwork/releases/download/v0.11.3/clockwork_v0.11.3.img':
        'ghcr.io/iqbal-lab-org/clockwork:latest' }"

    input:
    tuple val(meta), path(in_sam), path(ref_tsv)
    
    output:
    tuple val(meta), path("*{1,2}.fq.gz"), emit: reads
    tuple val(meta), path("*.counts.tsv"), emit: tsv
    
    script:
     def clean_now = task.ext.clean ? false : true
    """
     clockwork remove_contam $ref_tsv $in_sam ${meta.id}.decontam.counts.tsv ${meta.id}.decontam_1.fq.gz ${meta.id}.decontam_2.fq.gz
    
      ##  reshub addition, free the SSD disk ASAP
      if [[ ${params.clean_now} == true ]]; then
        rm \$(readlink ${in_sam})
      fi
    """
}
process REMOVE_CONTAM_MERGED{
    tag "$meta.id"
    label 'process_medium'
    
    memory { 5.GB * task.attempt }
	maxRetries 3
	errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/iqbal-lab-org/clockwork/releases/download/v0.11.3/clockwork_v0.11.3.img':
        'ghcr.io/iqbal-lab-org/clockwork:latest' }"

    input:
    tuple val(meta), path(reads), path(ref_tsv), val(ref_fa)
       

    output:
    tuple val(meta), path("*{1,2}.fq.gz"), emit: reads
    tuple val(meta), path("*.counts.tsv"), emit: tsv

    script:
     def threads = task.cpus
    """
     clockwork map_reads --threads $threads --unsorted_sam ${meta.id} $ref_fa ${meta.id}.sam ${reads[0]} ${reads[1]}
     clockwork remove_contam $ref_tsv ${meta.id}.sam ${meta.id}.decontam.counts.tsv ${meta.id}.decontam_1.fq.gz ${meta.id}.decontam_2.fq.gz
      
      ##  reshub addition, free the SSD disk ASAP
    if [[ ${params.clean_now} == true ]]; then
      rm -f ${meta.id}.sam
    fi

    """

}
process VARIANT_CALL {
    tag "$meta.id"
    label 'process_medium'

    memory { 5.GB * task.attempt }
	maxRetries 3
	errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
   
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/iqbal-lab-org/clockwork/releases/download/v0.11.3/clockwork_v0.11.3.img':
        'ghcr.io/iqbal-lab-org/clockwork:latest' }"

    input:
    tuple val(meta), path(reads), path(h37Rv_dir)
         

    output:
    tuple val(meta), path("*final.vcf"), emit: final_vcf
    tuple val(meta), path("*cortex.vcf"),optional:true, emit: cortex_vcf
    tuple val(meta), path("*samtools.vcf"),optional:true, emit: samtools_vcf

    publishDir { "${params.out_dir}/results/intermediate/${meta.id.tokenize('.')[0]}/cryptic"},
               mode: 'copy',
               pattern: '*.vcf'


    script:
     def clean_now = task.ext.clean ? false : true
    """
     clockwork variant_call_one_sample --sample_name ${meta.id} ./Ref.H37Rv var_call ${reads[0]} ${reads[1]}
     
     cp ./var_call/final.vcf ${meta.id}.final.vcf
     cp ./var_call/samtools.vcf ${meta.id}.samtools.vcf
     cp ./var_call/cortex.vcf ${meta.id}.cortex.vcf

      ##  reshub addition, free the SSD disk ASAP
     if [[ ${params.clean_now} == true ]]; then
        rm \$(readlink ${reads[0]})
        rm \$(readlink ${reads[1]})

        rm -f ./var_call/map.bam
        rm -f ./var_call/map.bai
     fi
    """
}
process PREDICT_DST{
    tag "$meta.id-$meta.cat "
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(catalog), path(ref_pkl) 
   
    output:
    tuple val(meta), path("*.effects.csv"), emit: effects
    tuple val(meta), path("*.mutations.csv"),optional:true, emit: mutations
    tuple val(meta), path("*.variants.csv"),optional:true, emit: variants
    tuple val(meta), path("*.json"),optional:true, emit: json

    publishDir { "${params.out_dir}/results/intermediate/${meta.id.tokenize('.')[0]}/cryptic"},
               mode: 'copy',
               pattern: '*.cryptic.*'
    script:
    
    """
    touch ${meta.id}.effects.csv 
    gnomonicus --vcf_file $vcf --genome_object $ref_pkl --catalogue_file $catalog  --output_dir . --json --progress
    mv *.effects.csv ${meta.id}.${meta.cat}.cryptic.effects.csv
    mv *.variants.csv ${meta.id}.${meta.cat}.cryptic.variants.csv
    mv *.mutations.csv ${meta.id}.${meta.cat}.cryptic.mutations.csv
    mv *.json ${meta.id}.${meta.cat}.cryptic.json
    """
} 