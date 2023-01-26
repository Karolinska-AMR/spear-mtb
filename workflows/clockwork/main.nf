 process MAP_READS {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

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
     
    """
     clockwork remove_contam $ref_tsv $in_sam ${meta.id}.decontam.counts.tsv ${meta.id}.decontam_1.fq.gz ${meta.id}.decontam_2.fq.gz
    """
}
process VARIANT_CALL {
    tag "$meta.id"
    label 'process_medium'
   
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/iqbal-lab-org/clockwork/releases/download/v0.11.3/clockwork_v0.11.3.img':
        'ghcr.io/iqbal-lab-org/clockwork:latest' }"

    input:
    tuple val(meta), path(reads)
    val h37Rv_dir
       

    output:
    tuple val(meta), path("*final.vcf"), emit: final_vcf
    tuple val(meta), path("*cortex.vcf"),optional:true, emit: cortex_vcf
    tuple val(meta), path("*samtools.vcf"),optional:true, emit: samtools_vcf

    script:
     
    """
     clockwork variant_call_one_sample --sample_name ${meta.id} $h37Rv_dir var_call ${reads[0]} ${reads[1]}
     cp ./var_call/final.vcf ${meta.id}.final.vcf
     cp ./var_call/samtools.vcf ${meta.id}.samtools.vcf
     cp ./var_call/cortex.vcf ${meta.id}.cortex.vcf
    """
}
process PREDICT_DST{
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(catalog), path(ref_pkl) 
    output:
    tuple val(meta), path("*.effects.csv"), emit: effects
    tuple val(meta), path("*.mutations.csv"),optional:true, emit: mutations
    tuple val(meta), path("*.variants.csv"),optional:true, emit: variants

    script:
    
    """
    python ${baseDir}/bin/sp3predict.py --prefix ${meta.id} --vcf_file $vcf --catalogue_file $catalog --genome_object $ref_pkl --progress --ignore_vcf_status --ignore_vcf_filter --debug

    """
} 