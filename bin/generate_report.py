from glob import glob
import os
import json
import argparse
from datetime import datetime


def convert_drugname(abbr):

    drug_mapper = {"AMC": "AMOXICILIN-CLAVULANATE", "AMI": "AMIKACIN", "AMX": "AMOXICILIN", "AZM": "AZITHROMYCIN", "BDQ": "BEDAQUILINE", "CAP": "CAPREOMYCIN",
                   "CFZ": "CLOFAZIMINE", "CIP": "CIPROFLOXACIN", "CLR": "CLARITHROMYCIN", "CYC": "CYCLOSERINE", "DCS": "D-CYCLOSERINE", "DLM": "DELAMANID",
                   "EMB": "ETHAMBUTOL", "ETH": "ETHIONAMIDE", "ETP": "ERTAPENEM", "FQS": "FLUOROQUINOLONES", "GEN": "GENTAMICIN", "GFX": "GATIFLOXACIN",
                   "IMI": "IMIPENEM", "INH": "ISONIAZID", "KAN": "KANAMYCIN", "LEV": "LEVOFLOXACIN", "LZD": "LINEZOLID", "MEF": "MEFLOQUINE", "MPM": "MEROPENEM",
                   "MXF": "MOXIFLOXACIN", "OFX": "OFLOXACIN", "PAN": "PRETOMANID", "PAS": "PARA-AMINOSALICYLIC_ACID", "PTO": "PROTHIONAMIDE", "PZA": "PYRAZINAMIDE", "RFB": "RIFABUTIN",
                   "RIF": "RIFAMPICIN", "STM": "STREPTOMYCIN", "STX": "SITAFLOXACIN", "SXT": "COTRIMOXAZOLE", "SZD": "SUTEZOLID", "TRD": "TERIZIDONE", "TZE": "THIOACETAZONE"}
    try:
        return drug_mapper[abbr]
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

    mapper = {
        "fgd1": [
            ["delamanid", "R"],
            ["clofazimine", "U"]
        ],
        "ethA": [
            ["ethionamide", "R"]
        ],
        "katG": [
            ["isoniazid", "R"]
        ],
        "pncA": [
            ["pyrazinamide", "R"]
        ],
        "bacA": [
            ["amikacin", "U"],
            ["capreomycin", "U"],
            ["kanamycin", "U"],
            ["streptomycin", "U"]
        ],
        "eis": [
            ["amikacin", "U"],
            ["kanamycin", "U"]
        ],
        "whiB6": [
            ["amikacin", "U"],
            ["capreomycin", "U"],
            ["kanamycin", "U"]
        ],
        "Rv1979c": [
            ["bedaquiline", "U"],
            ["clofazimine", "U"]
        ],
        "ndh": [
            ["delamanid", "U"],
            ["ethionamide", "U"],
            ["isoniazid", "U"]
        ],
        "embR": [
            ["ethambutol", "U"]
        ],
        "glpK": [
            ["ethambutol", "U"],
            ["isoniazid", "U"],
            ["levofloxacin", "U"],
            ["moxifloxacin", "U"],
            ["rifampicin", "U"],
            ["streptomycin", "U"]
        ],
        "Rv2752c": [
            ["ethambutol", "U"],
            ["isoniazid", "U"],
            ["levofloxacin", "U"],
            ["moxifloxacin", "U"],
            ["rifampicin", "U"]
        ],
        "mshA": [
            ["ethionamide", "U"],
            ["isoniazid", "U"]
        ],
        "Rv3083": [
            ["ethionamide", "U"]
        ],
        "ahpC": [
            ["isoniazid", "U"]
        ],
        "Rv1258c": [
            ["isoniazid", "U"],
            ["pyrazinamide", "U"],
            ["streptomycin", "U"]
        ],
        "PPE35": [
            ["pyrazinamide", "U"]
        ],
        "Rv3236c": [
            ["pyrazinamide", "U"]
        ]
    }

    try:
        return mapper[gene]
    except KeyError:
        return []


def get_base_dict(sample_id):

    return {
        "seqid": sample_id,
        "category": None,
        'catalogs':[],
        "pct_mapped_reads": None,
        "num_mapped_reads": None,
        "region_median_depth": None,
        "genome_median_depth": None,
        "lineage": {},
        "variants": [],
        "unverified_mutations": {},
        "del_genes": [],
        "resist_drug": {}
    }


def get_catalog_name(cat_name):
        
        if cat_name.startswith('who'):
            catalog = "WHO-UCN-2023.5"
        elif cat_name.startswith('tbdb'):
            catalog = "TB-Profiler"
        elif cat_name.lower().startswith('cryptic'):
            catalog = 'CRyPTIC'
        else:
            catalog = cat_name
        
        return catalog


def collect_tbprofilers(rep_dir, report_dict):

    # Extract tbprofiler dictionary
    tbp_out_lst = sorted(
        glob(os.path.join(rep_dir, "*.tbprofiler.results.json")))

    for fj in tbp_out_lst:
     
        sid = os.path.basename(fj).split('.')[0]
        cat = os.path.basename(fj).split('.')[1]

        catalog = get_catalog_name(cat)

        if sid not in report_dict:
            report_dict[sid] = get_base_dict(sid)

        ptr = report_dict[sid]
        ptr['catalogs'].append(catalog)
        with open(fj) as hdl:
            report_json = json.load(hdl)
            if "qc" in report_json:
                ptr_qc = report_json["qc"]
                ptr["pct_mapped_reads"] = ptr_qc["pct_reads_mapped"]
                ptr["num_mapped_reads"] = ptr_qc["num_reads_mapped"]
                ptr["region_median_depth"] = ptr_qc["region_median_depth"]
                ptr["genome_median_depth"] = ptr_qc["genome_median_depth"]

            if "main_lin" in report_json:

                tmp_ptr = {"lin": [], "family": [],
                           "spoligotype": [], 
                            "rd": [], "frac": [] }

                prv_lin = ""
                lineages = sorted(
                    report_json["lineage"], key=lambda x: x['lin'])
                for lng in lineages[::-1]:
                    if not prv_lin.startswith(lng['lin']):
                        tmp_ptr['frac'].append(round(lng['frac'], 2))
                        tmp_ptr['lin'].append(lng['lin'])
                        tmp_ptr['family'].append(lng['family'])
                        tmp_ptr['spoligotype'].append(lng['spoligotype'])
                        tmp_ptr['rd'].append(lng['rd'])
                    prv_lin = lng['lin']

                is_mixed = len(tmp_ptr["lin"]) > 1
                for k, v in tmp_ptr.items():
                    tmp_ptr[k] = '|'.join(map(str, v)).replace('None', '')
                tmp_ptr["is_mixed"] = is_mixed
                ptr['lineage'] = tmp_ptr

            ptr_var = ptr['variants']
            
            for rvr in report_json["dr_variants"]:
                for drg in rvr["drugs"]:
                    name = str(drg["drug"]).upper()
                    if drg["confers"] != "sensitive":
                        if drg["confers"]== 'resistance':
                            pred = "R"
                        elif drg["confers"]== 'uncertain':
                            pred = "U"
                        else:
                            pred = '?'

                        ptr_var.append({'catalog': catalog, "gene": rvr['gene'],'drug': name,
                               "variant": rvr['change'], "pred": pred,
                               'freq': round(rvr['freq'], 2), 'depth': rvr['depth']})

                        if pred == 'R':
                            if name not in ptr['resist_drug']:
                                ptr['resist_drug'][name] = set()
                            ptr['resist_drug'][name].add(catalog)

            for rvr in report_json["other_variants"]:
                
                if "annotation" not in rvr:
                    continue

                for annot in rvr["annotation"]:
                    name = str(annot["drug"]).upper()
                    if "who_confidence" in annot:
                        grade = annot["who_confidence"]
                    elif "grade" in annot:
                        grade = annot["grade"]
                    else:
                        print(f'UNKNOWN FORMAT: {annot}')
                        continue
                    if grade.lower().find('uncertain')>-1:    
                        ptr_var.append({'catalog': catalog, "gene": rvr['gene'],'drug': name,
                                "variant": rvr['change'], "pred": 'U',
                                'freq': round(rvr['freq'], 2), 'depth': rvr['depth']})

            del_genes = get_gene_del(report_json["qc"]["region_qc"])
            ptr_dg = ptr["del_genes"]
            for gene, drugs in del_genes.items():
                for drg,pred in drugs:
                    ptr_dg.append({'catalog': catalog, "gene": gene,'pred': pred,'drug':drg,
                                  "variant": "feature_ablation", 'freq': None,'depth': None})

            ptr['unverified_mutations'] = get_unverfied_regions(ptr['unverified_mutations'],
                                                                report_json["qc"]["missing_positions"],
                                                                del_genes, catalog)

    return report_dict


def get_gene_del(region_qc):
    del_genes = {}
    for item in region_qc:
        reg = item['region']
        gene = convert_locusname(reg)
        if item['median_depth'] <= 5:
            del_genes[gene] = gene_del_drugs(gene)

    return del_genes


def get_unverfied_regions(unverified_variants, missing_pos, missing_genes, catalog):
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
            # it's deletion, not missing R
            if var.find('del') > -1:
                continue
            key = f"{gene}_{var}"
            if key not in unverified_variants:
                unverified_variants[key] = (
                    {'catalog': catalog, "gene": gene, "locus": locus, "variant": var})
    return unverified_variants


def collect_cryptics(rep_dir, report_dict):
    cryp_out_lst = glob(os.path.join(rep_dir, "*cryptic.json"))

    for fj in cryp_out_lst:
        sid = os.path.basename(fj).split('.')[0]
        cat = os.path.basename(fj).split('.')[1]
        
        catalog = get_catalog_name(cat)
        
        if sid not in report_dict:
            report_dict[sid] = get_base_dict(sid)

        ptr = report_dict[sid]
        ptr['catalogs'].append(catalog)
        with open(fj) as hdl:
            res_dict = json.load(hdl)
            res_dict = res_dict["data"]

        ptr_var = ptr['variants']
        if "EFFECTS" in res_dict:
            for drg, mut_list in res_dict["EFFECTS"].items():
                drg = convert_drugname(drg)

                for mut in mut_list:

                    if "PHENOTYPE" in mut:
                        continue

                    if mut['PREDICTION'] != "S":
                        ptr_var.append({'catalog': catalog, "gene": mut['GENE'],
                                        'pred':mut['PREDICTION'],'drug':drg,
                                        "variant": mut['MUTATION'],
                                        'freq': None,'depth':None})

                        if mut['PREDICTION'] == "R":
                            if drg not in ptr['resist_drug']:
                                ptr['resist_drug'][drg] = set()
                            ptr['resist_drug'][drg].add(catalog)

    return report_dict


def classify_tb_resistance(resistant_drugs):

    if len(resistant_drugs) == 0:
        return 'DS'

    # Define the key drugs for various TB resistance categories
    first_line_drugs = {'ISONIAZID', 'RIFAMPICIN',
                        'ETHAMBUTOL', 'PYRAZINAMIDE'}
    second_line_injectables = {'AMIKACIN', 'KANAMYCIN', 'CAPREOMYCIN'}
    fluoroquinolones = {'FLUOROQUINOLONE', 'MOXIFLOXACIN', 'OFLOXACIN',
                        'LEVOFLOXACIN', 'GATIFLOXACIN', 'CIPROFLOXACIN', 'SITAFLOXACIN'}
    group_A = {'MOXIFLOXACIN', 'LEVOFLOXACIN', 'BEDAQUILINE', 'LINEZOLID'}
    # Check for various resistance types
    is_INH_resistant = 'ISONIAZID' in resistant_drugs
    is_RIF_resistant = 'RIFAMPICIN' in resistant_drugs
    is_first_line_resistant = first_line_drugs.intersection(resistant_drugs)
    is_second_line_injectable_resistant = second_line_injectables.intersection(
        resistant_drugs)
    is_fluoroquinolone_resistant = fluoroquinolones.intersection(
        resistant_drugs)
    is_group_a_resistant = group_A.intersection(resistant_drugs)

    if len(is_first_line_resistant) == 1:
        # Resistance to a single first-line anti-TB drug
        return 'Mono'
    elif is_INH_resistant and is_RIF_resistant:
        if is_fluoroquinolone_resistant:
            if is_group_a_resistant:
                return 'XDR'
            return 'Pre-XDR'
        return 'MDR'
    # elif is_INH_resistant and not is_RIF_resistant:
    #     return 'IR'
    # elif is_RIF_resistant and not is_INH_resistant:
    #     return 'RR'
    elif len(is_first_line_resistant) >= 2:
        return 'Poly'
    # elif is_INH_resistant or is_RIF_resistant:
    #     return 'Pre-MDR'
    else:
        return 'other'


def generate(in_dir, prefix):
    report_dict = {}
    collect_tbprofilers(in_dir, report_dict)
    collect_cryptics(in_dir, report_dict)

    for sid, content in report_dict.items():
        report_dict[sid]['category'] = classify_tb_resistance(
            list(content['resist_drug'].keys()))
        report_dict[sid]['resist_drug'] = list(content['resist_drug'])

    djson = json.dumps(report_dict)

    with open(f"{prefix}.json", 'w') as hdl:
        hdl.write(djson)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_dir", required=True,
                        help="Iput directory containing results of ")

    parser.add_argument("--prefix", required=True,
                        help="Ticket number")

    options = parser.parse_args()

    generate(options.in_dir, options.prefix)
