nextflow.enable.dsl=2

// Include workflows and modules


include { REMOVE_CONTAM_MERGED} from "$baseDir/workflows/clockwork/main"
include { VARIANT_CALL}  from "$baseDir/workflows/clockwork/main"
include { PREDICT_DST}   from "$baseDir/workflows/clockwork/main"
include { MTB_FINDER }           from "$baseDir/workflows/myco_miner/main"

include { TBPROFILER_PROFILE }     from "$baseDir/workflows/tbprofiler/profile/main"
include { TBPROFILER_PROFILE_EXTERNAL } from "$baseDir/modules/tbprofiler_external/main"
include { GENERATE_REPORT }        from "$baseDir/modules/utility/main"
include { ARCHIVE_RAW_SEQ }       from "$baseDir/modules/utility/main"

// PARAMETERS
params.input_dir        = ""
params.out_dir          = ""
params.k2_db            = "k2_myco"
params.assets_dir       = "${baseDir}/assets"
params.archive_input    = params.archive_input ?: false
params.archive_dir      =params.archive_dir ?: ""
params.ticket = params.ticket ?: "ticket"
params.skip_cryptic = params.skip_cryptic in [true, 'true', ''] ? true : false

def h37Rv_dir     = "${params.assets_dir}/Ref.H37Rv"
def kraken2_db    = "${params.assets_dir}/kraken2/${params.k2_db}"
def template_report = "${params.assets_dir}/report/report-template.html"
workflow {

    // Ensure output directory exists
    file(params.out_dir).mkdir()

    // Check archive parameters
   if (params.archive_input && !params.archive_dir) {
      error "archive_dir must be set if archive_input is true"
   }

    //  READ INPUTS 
    Channel
        .fromFilePairs("${params.input_dir}/*_{1,2}.{fq.gz,fastq.gz}")
        .concat(Channel.fromFilePairs("${params.input_dir}/*/*_{1,2}.{fq.gz,fastq.gz}"))
        .set { input_reads_ch }
      

   input_reads_ch
        .count()
        .map{count ->
            if (count==0) {
                error "No input data provided. Please check the input directory: ${params.input_dir}"
            }
        }
    
    //  SPECIES FILTERING: Making sure we only process MTB genomes 
    MTB_FINDER(input_reads_ch, kraken2_db)
    MTB_FINDER.out.mtbs_ch.set { mtb_reads_ch }

    //  TBProfiler 
    mtb_reads_ch.map { [ [ id: "${it[0]}.tbdb.tbprofiler", single_end: false ], it[1] ] }
                .set { tbprofiler_input_ch }

    TBPROFILER_PROFILE(tbprofiler_input_ch)

    // TBProfiler External Catalogue
    TBPROFILER_PROFILE.out.bam.map { [ [ id: "${it[1].simpleName}.${params.catalogue_version}.tbprofiler" ], it[1] ] }
                .combine(Channel.fromPath("${params.assets_dir}/catalogues/NC_000962.3/who2023v07"))
                .set { tbprofiler_external_ch }

    TBPROFILER_PROFILE_EXTERNAL(tbprofiler_external_ch)

    //  CRYPTIC WORKFLOW 
   if (!params.skip_cryptic) {
      Channel.fromPath("${params.assets_dir}/Ref.remove_contam/*.tsv")
            .set { contam_tsv_ch }

      Channel.fromPath("${params.assets_dir}/Ref.remove_contam/*.fa")
            .set { contam_fa_ch }

      mtb_reads_ch
         .combine(contam_tsv_ch)
         .combine(contam_fa_ch)
         .map { [ [ id: "${it[0]}.cryptic" ], it[1], it[2], it[3] ] }
         .set { cryptic_input_ch }

      REMOVE_CONTAM_MERGED(cryptic_input_ch)

      //  VARIANT CALLING 
      REMOVE_CONTAM_MERGED.out.reads
            .combine(Channel.fromPath(h37Rv_dir))
            .set { variant_call_input_ch }

      VARIANT_CALL(variant_call_input_ch)

      //  DST PREDICTION 
      Channel.fromPath("${params.assets_dir}/catalogues/*/*.csv")
            .set { catalog_csv_ch }

      Channel.fromPath("${params.assets_dir}/catalogues/*/*.gz")
            .set { refpkl_gz_ch }

      VARIANT_CALL.out.final_vcf
         .combine(catalog_csv_ch)
         .combine(refpkl_gz_ch)
         .map { [ [ id: it[1].simpleName, cat: it[2].simpleName ], it[1], it[2], it[3] ] }
         .set { dst_input_ch }

      PREDICT_DST(dst_input_ch)

      //  COLLECT OUTPUT JSONs 
      PREDICT_DST.out.json
         .concat(TBPROFILER_PROFILE.out.json)
         .concat(TBPROFILER_PROFILE_EXTERNAL.out.json)
         .map { it[1] }
         .flatten()
         .set { all_jsons_ch }
   } else {
      // If not running cryptic, just collect TBProfiler outputs
      TBPROFILER_PROFILE.out.json
         .concat(TBPROFILER_PROFILE_EXTERNAL.out.json)
         .map { it[1] }
         .flatten()
         .set { all_jsons_ch }
   }

    //  ARCHIVE IF NEEDED 
    if (params.archive_input) {
        mtb_reads_ch.map { it[0] }
                    .join(PREDICT_DST.out.json.concat(TBPROFILER_PROFILE.out.json).map { it[0].id.split("\\.")[0] })
                    .set { archived_reads_ch }

        ARCHIVE_RAW_SEQ(archived_reads_ch, params.archive_dir, 'mv')
    }

    //  GENERATE REPORT 
    GENERATE_REPORT(all_jsons_ch.collect(), params.ticket, template_report)

}

//  OPTIONAL: Finalize 
workflow.onComplete {
   println "Pipeline completed. Output available at: ${params.out_dir}"
   println "Ticket: ${params.ticket}"
}
