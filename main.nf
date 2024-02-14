nextflow.enable.dsl=2

include {MAP_READS as mpr} from "$baseDir/workflows/clockwork/main"
include {REMOVE_CONTAM as rmc} from "$baseDir/workflows/clockwork/main"
include {REMOVE_CONTAM_MERGED as rmc_mrg} from "$baseDir/workflows/clockwork/main"
include {VARIANT_CALL as vrc} from "$baseDir/workflows/clockwork/main"
include {PREDICT_DST as prd} from "$baseDir/workflows/clockwork/main"
include {MTB_FINDER} from "$baseDir/workflows/myco_miner/main"

include {TBPROFILER_PROFILE as tbp} from "$baseDir/workflows/tbprofiler/profile/main"
include {TBPROFILER_PROFILE_EXTERNAL as tbp_ext} from "$baseDir/modules/tbprofiler_external/main"
include {GENERATE_REPORT as grp} from "$baseDir/modules/utility/main"
include {ARCHIVE_RAW_SEQ as arch} from "$baseDir/modules/utility/main"

params.input_dir = ""

def assets_dir = params.assets_dir ?:"${baseDir}/assets"
def h37Rv_dir = "${assets_dir}/Ref.H37Rv";
def out_dir = params.out_dir
def k2_db='k2_myco'

workflow{


   file(out_dir).mkdir()

   reads_ch = Channel.fromFilePairs("${params.input_dir}/*_{1,2}.fastq.gz")
                     .concat(Channel.fromFilePairs("${params.input_dir}/*/*_{1,2}.fastq.gz"));

   ref_fa = Channel.fromPath("${assets_dir}/Ref.remove_contam/*.fa");
   
   //find species and continue with MTBs
   
    kraken2_db = "${assets_dir}/kraken2/${k2_db}"
    MTB_FINDER(reads_ch,kraken2_db)
    input_ch = MTB_FINDER.out.mtbs_ch
  
   // Running TBPROFILER
   
   tbpr_ch = input_ch.map{it->[[id:"${it[0]}.tbdb.tbprofiler",single_end:false],it[1]]}
   tbp(tbpr_ch,out_dir)
  
   bam_ch = tbp.out.bam.map{it->[[id:"${it[1].simpleName}.who2023v5.tbprofiler"],it[1]]}
                .combine(Channel.fromPath("${assets_dir}/catalogues/NC_000962.3/WHO-2023.5"))
   tbp_ext(bam_ch,out_dir)

   // Running CRyPTIC workflow
   contam_ref_ch = Channel.fromPath("${assets_dir}/Ref.remove_contam/*.tsv");

//    cryptic_ch =  input_ch.combine(ref_fa).map{it->[[id:"${it[0]}.cryptic"],it[1],it[2]]}  
//    mpr(cryptic_ch) 
//    rmc(mpr.out.sam.combine(contam_ref_ch))

   cryptic_ch =  input_ch.combine(contam_ref_ch).combine(ref_fa).map{it->[[id:"${it[0]}.cryptic"],it[1],it[2]]}
   rmc_mrg(cryptic_ch)

   vrc(rmc_mrg.out.reads.combine(Channel.fromPath(h37Rv_dir)),out_dir)

   catalog_ch = Channel.fromPath("${assets_dir}/catalogues/*/*.csv")
   refpkl_ch = Channel.fromPath("${assets_dir}/catalogues/*/*.gz")

   ch = vrc.out.final_vcf.combine(catalog_ch).combine(refpkl_ch)
           .map{it->[[id:it[1].simpleName,cat:it[2].simpleName],it[1],it[2],it[3]]}

   prd(ch,out_dir)
    
   ser_ch = prd.out.json
          .concat(tbp.out.json)
          .concat(tbp_ext.out.json).map{it->it[1]}.flatten()
          
    // Archive the successfully executed Illumina PE-reads
    arc_ch = input_ch.join(prd.out.json.concat(tbp.out.json).map{it-> it[0].id.split("\\.")[0]})

    arch(arc_ch,params.archive,'mv')

    grp(ser_ch.collect(),params.ticket)
    grp.out.json.collectFile(storeDir:"${out_dir}/reports")
} 
