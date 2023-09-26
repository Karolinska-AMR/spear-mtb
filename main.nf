nextflow.enable.dsl=2

include {MAP_READS as mpr} from "$baseDir/workflows/clockwork/main"
include {REMOVE_CONTAM as rmc} from "$baseDir/workflows/clockwork/main"
include {VARIANT_CALL as vrc} from "$baseDir/workflows/clockwork/main"
include {PREDICT_DST as prd} from "$baseDir/workflows/clockwork/main"
include {MTB_FINDER} from "$baseDir/workflows/myco_miner/main"

include {TBPROFILER_PROFILE as tbp} from "$baseDir/workflows/tbprofiler/profile/main"
include {GENERATE_REPORT as grp} from "$baseDir/modules/utility/main"


params.input_dir = ""

def assets_dir = "${baseDir}/assets"
def h37Rv_dir = "${assets_dir}/Ref.H37Rv";
def out_dir = params.out_dir ?: "${baseDir}/out"
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
    input_ch.view()
   // Running TBPROFILER
   tbpr_ch = input_ch.map{it->[[id:"${it[0]}.tbprofiler",single_end:false],it[1]]}
   tbp(tbpr_ch)

   // Running CRyPTIC workflow
   contam_ref_ch = Channel.fromPath("${assets_dir}/Ref.remove_contam/*.tsv");
   cryptic_ch =  input_ch.combine(ref_fa).map{it->[[id:"${it[0]}.cryptic"],it[1],it[2]]}
   
   mpr(cryptic_ch) 
   rmc(mpr.out.sam.combine(contam_ref_ch))
   vrc(rmc.out.reads,h37Rv_dir)
   
   
   vrc.out.final_vcf.map{it-> it[1]}.flatten().collectFile(storeDir:out_dir)

   catalog_ch = Channel.fromPath("${assets_dir}/catalogues/*/*.csv")
   refpkl_ch = Channel.fromPath("${assets_dir}/catalogues/*/*.gz")

   ch = vrc.out.final_vcf.combine(catalog_ch).combine(refpkl_ch)
           .map{it->[[id:it[1].simpleName,cat:it[2].simpleName],it[1],it[2],it[3]]}
   prd(ch)
    
    
   ser_ch = prd.out.effects
          .concat(prd.out.mutations)
          .concat(prd.out.variants)
          .concat(prd.out.json)
          .concat(tbp.out.txt)
          .concat(tbp.out.json)
          .map{it->it[1]}
          .concat(MTB_FINDER.out.species_tsvs)
          .flatten()
          .collectFile(storeDir:out_dir)

    grp(ser_ch.collect(),"${assets_dir}/report/report-template.html")
    grp.out.html.collectFile(storeDir:out_dir)
} 

workflow generate_report{
    json_ch = Channel.fromPath("${out_dir}/*.json").collect()
    grp(json_ch,"${assets_dir}/report/report-template.html")
    grp.out.html.collectFile(storeDir:out_dir)
}