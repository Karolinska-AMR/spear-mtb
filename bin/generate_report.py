from glob import glob
import os
import json
import argparse
from datetime import datetime


def convert_drugname(abbr):

    drug_mapper = {"AMC": "AMOXICILIN-CLAVULANATE", "AMI": "AMIKACIN", "AMX": "AMOXICILIN", "AZM": "AZITHROMYCIN", "BDQ": "BEDAQUILINE", "CAP": "CAPREOMYCIN",
                   "CFZ": "CLOFAZIMINE", "CIP": "CIPROFLOXACIN", "CLR": "CLARITHROMYCIN", "CYC": "CYCLOSERINE", "DCS": "D-CYCLOSERINE", "DLM": "DELAMANID",
                   "EMB": "ETHAMBUTOL", "ETH": "ETHIONAMIDE", "ETP": "ERTAPENEM", "FQS": "FLUOROQUINOLONE", "GEN": "GENTAMICIN", "GFX": "GATIFLOXACIN",
                   "IMI": "IMIPENEM", "INH": "ISONIAZID", "KAN": "KANAMYCIN", "LEV": "LEVOFLOXACIN", "LZD": "LINEZOLID", "MEF": "MEFLOQUINE", "MPM": "MEROPENEM",
                   "MXF": "MOXIFLOXACIN", "OFX": "OFLOXACIN", "PAN": "PRETOMANID", "PAS": "PAS", "PTO": "PROTHIONAMIDE", "PZA": "PYRAZINAMIDE", "RFB": "RIFABUTIN",
                   "RIF": "RIFAMPICIN", "STM": "STREPTOMYCIN", "STX": "SITAFLOXACIN", "SXT": "COTRIMOXAZOLE", "SZD": "SUTEZOLID", "TRD": "TERIZIDONE", "TZE": "THIOACETAZONE"}
    try:
        return drug_mapper[abbr].lower()
    except KeyError:
        return abbr


def convert_locusname(name):
    locus_mapper = {"Rv0005": "gyrB", "Rv0006": "gyrA", "Rv0407": "fgd1", "Rv0486": "mshA", "Rv0667": "rpoB", "Rv0668": "rpoC", "Rv0678": "mmpR5", "Rv0682": "rpsL", "Rv0701": "rplC",
                    "Rv1173": "fbiC", "Rv1267c": "embR", "Rv1305": "atpE", "rrs": "rrs", "rrl": "rrl", "Rv1483": "fabG1", "Rv1484": "inhA", "Rv1630": "rpsA", "Rv1694": "tlyA", "Rv1908c": "katG", "Rv2043c": "pncA",
                    "Rv2245": "kasA", "Rv2416c": "eis", "Rv2428": "ahpC", "Rv2447c": "folC", "Rv2535c": "pepQ", "Rv2671": "ribD", "Rv2754c": "thyX", "Rv2764c": "thyA", "Rv2780": "ald", "Rv2983": "fbiD", "Rv3261": "fbiA",
                    "Rv3262": "fbiB", "Rv3423c": "alr", "Rv3547": "ddn", "Rv3601c": "panD", "Rv3793": "embC", "Rv3794": "embA", "Rv3795": "embB", "Rv3806c": "ubiA", "Rv3854c": "ethA", "Rv3855": "ethR", "Rv3919c": "gid"}

    try:
        return locus_mapper[name]
    except KeyError:
        return name


def gene_del_drugs(gene):

    mapper = {"katG": ["isoniazid"],
              "pncA": ["pyrazinamide"],
              "fgd1": ["delamanid"],
              "ethA": ["ethionamide"]
              }
    try:
        return mapper[gene]
    except KeyError:
        return []


def collect_tbprofilers(rep_dir, report_dict):

    # Extract tbprofiler dictionary
    tbp_out_lst = sorted(
                 glob(os.path.join(rep_dir, "*.tbprofiler.results.json")))
    
    for fj in tbp_out_lst:
        sid = os.path.basename(fj).split('.')[0]
        cat = os.path.basename(fj).split('.')[1]
       
        if cat.startswith('who'):
            catalog = "WHO-UCN-2023.5"
        elif cat.startswith('tbdb'):
            catalog = "TB-Profiler"
        else:
            catalog="unknown"

        if sid not in report_dict:
            report_dict[sid] = {"seqid": sid, "lineage": {},
                                "variants": [], "unverified_mutations": {},
                                "del_genes": []}
            
        ptr = report_dict[sid]
        with open(fj) as hdl:
            report_json = json.load(hdl)
            # Directly extracting the sublineage
            if "main_lin" in report_json:
                pp = report_json["main_lin"].split(';')
                if len(pp) == 1:
                    if len(report_json["lineage"]):
                        # selecting the last sublineage
                        ptr["lineage"] = report_json["lineage"][-1]
                        ptr["lineage"]["is_mixed"] = False
                elif len(pp) > 1:
                    ptr["lineage"] = {"lin": report_json["main_lin"], "family": None,
                                      "spoligotype": None, "rd": None, "frac": None}
                    ptr["lineage"]["is_mixed"] = True
                  

            ptr_var = ptr['variants']
            variations = report_json["dr_variants"] +  report_json["other_variants"]
            for rvr in variations:   
                ptr_var.append({'catalog':catalog, "gene": rvr['gene'],"variant":rvr['change'],'freq':round(rvr['freq'], 2)})
           
           
            del_genes = get_gene_del(report_json["qc"]["region_qc"])
            ptr_dg = ptr["del_genes"]
            for gene, drugs in del_genes.items():
                for drg in drugs:
                    if drg not in ptr_dg:
                        ptr_dg[drg] = []
                    ptr_dg.append({'catalog':catalog, "gene": gene,"variant":"feature_ablation",'freq':None})
            
            
            ptr['unverified_mutations'] = get_unverfied_regions(ptr['unverified_mutations'],
                                  report_json["qc"]["missing_positions"],
                                  del_genes,catalog)

    return report_dict


def get_gene_del(region_qc):
    del_genes = {}
    for item in region_qc:
        reg = item['region']
        gene = convert_locusname(reg)
        if item['median_depth'] <= 5:
            del_genes[gene] = gene_del_drugs(gene)

    return del_genes


def get_unverfied_regions( unverified_variants, missing_pos, missing_genes,catalog):
    ctrl_rec = set()
    for pos in missing_pos:
        gene = pos['gene']
        locus = pos['locus_tag']
        # drugs = pos['drugs'].split(',')
        variants = pos['variants'].split(',')
        is_deleted = False
        if gene in missing_genes:
            is_deleted = True
            variants = [f"feature_ablation"]

        # for drg in drugs:
        if is_deleted:
            key = f"{gene}"
            if key in ctrl_rec:
                continue
            ctrl_rec.add(key)

        for var in variants:
            #it's deletion, not missing R
            if var.find('del')>-1:
                continue
            key = f"{gene}_{var}"
            if key not in unverified_variants:
                unverified_variants[key] = (
                    {'catalog':catalog,'position':pos['position'], "gene": gene, "locus": locus, "variant": var})
    return unverified_variants


def collect_cryptics(rep_dir, report_dict):
    cryp_out_lst = glob(os.path.join(rep_dir, "*cryptic.json"))

    for fj in cryp_out_lst:
        sid = os.path.basename(fj).split('.')[0]
        catalog = os.path.basename(fj).split('.')[1]

        if sid not in report_dict:
            report_dict[sid] = {"seqid": sid, "lineage": {},
                                "variants": [], "unverified_mutations": {},
                                "del_genes": []}
            
        ptr = report_dict[sid]
        with open(fj) as hdl:
            res_dict = json.load(hdl)
            res_dict = res_dict["data"]

        # if catalog not in ptr["res_var"]:
        #     ptr["res_var"][catalog] = {}
        # ptr_cat = ptr["res_var"][catalog]
        ptr_var = ptr['variants']
        if "EFFECTS" in res_dict:
            for drg, mut_list in res_dict["EFFECTS"].items():
                drg = convert_drugname(drg)

                for mut in mut_list:

                    if "PHENOTYPE" in mut:
                        continue

                    if mut['PREDICTION'] != "S":               
                        ptr_var.append({'catalog':catalog, "gene": mut['GENE'],
                                        "variant":f"{mut['GENE']}@{mut['MUTATION']}",
                                                'freq':None})
                            

    return report_dict


def generate(in_dir, template):
    report_dict = {}
    collect_tbprofilers(in_dir, report_dict)
    collect_cryptics(in_dir, report_dict)
    djson = json.dumps(report_dict)  

    with open(f"spear-mtb_report_summary.json",'w') as hdl:
        hdl.write(djson)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_dir", required=True,
                        help="Iput directory containing results of ")
    parser.add_argument("--html", required=True,
                        help="The report html template")

    options = parser.parse_args()

    generate(options.in_dir, options.html)
