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


def convert_gene_locus(name, reverse=False):
    data = {"gyrB": "Rv0005", "gyrA": "Rv0006", "fgd1": "Rv0407", "mshA": "Rv0486", "rpoB": "Rv0667",
                    "rpoC": "Rv0668", "mmpR5": "Rv0678", "rpsL": "Rv0682", "rplC": "Rv0701", "fbiC": "Rv1173",
                    "embR": "Rv1267c", "atpE": "Rv1305", "rrs": "rrs", "rrl": "rrl", "fabG1": "Rv1483",
                    "inhA": "Rv1484", "rpsA": "Rv1630", "tlyA": "Rv1694", "katG": "Rv1908c", "pncA": "Rv2043c",
                    "kasA": "Rv2245", "eis": "Rv2416c", "ahpC": "Rv2428", "folC": "Rv2447c", "pepQ": "Rv2535c",
                    "ribD": "Rv2671", "thyX": "Rv2754c", "thyA": "Rv2764c", "ald": "Rv2780", "fbiD": "Rv2983",
                    "fbiA": "Rv3261", "fbiB": "Rv3262", "alr": "Rv3423c", "ddn": "Rv3547", "panD": "Rv3601c",
                    "embC": "Rv3793", "embA": "Rv3794", "embB": "Rv3795", "ubiA": "Rv3806c", "ethA": "Rv3854c",
                    "ethR": "Rv3855", "gid": "Rv3919c"
            }
    if reverse:
        data = {v: k for k, v in data.items()}

    return data.get(name, name)


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
    return mapper.get(gene, [])


def get_base_dict(sample_id):

    return {
        "seqid": sample_id,
        "category_who": None,
        'analyzed_catalogs': [],
        "percent_reads_mapped": None,
        "num_mapped_reads": None,
        "target_median_depth": None,
        "genome_median_depth": None,
        "lineage": {},
        "variants": [],
        "unverified_mutations": {},
        "del_genes": [],
        "resist_drug": [],
        "notes": [],
    }


def get_catalog_repo():
    return {'who':"WHO__v2023.07",'tbdb':"TB-Profiler__v6.4.0",'cryptic':"CRyPTIC__v1-311-6627d81_WHOv2-7577f14"}

def get_catalog_fullname(name):
    repo = get_catalog_repo()
    return repo.get(name,name)

def collect_tbprofilers(rep_dir, report_dict):

    # Extract tbprofiler dictionary
    tbp_out_lst = sorted(
        glob(os.path.join(rep_dir, "*.tbprofiler.results.json")))

    for fj in tbp_out_lst:

        sid = os.path.basename(fj).split('.')[0]
        cat = os.path.basename(fj).split('.')[1]
        catch_resis_drg = False
        if cat.startswith('who'):
            cat = 'who'
            catch_resis_drg = True

        catalog = get_catalog_fullname(cat)

        if sid not in report_dict:
            report_dict[sid] = get_base_dict(sid)

        ptr = report_dict[sid]
        
        with open(fj) as hdl:
            report_json = json.load(hdl)

            # # directly from the json file
            # if "pipeline" in report_json:
            #     catalog = report_json["pipeline"]["db_version"]["name"]
            #     if catalog.startswith('who'):
            #         catalog = 'WHO_v2023.7'
            #     elif catalog.startswith('tbdb'):
            #         catalog = 'TB-Profiler_v'+ report_json["pipeline"]['software_version']
            
            ptr['analyzed_catalogs'].append(catalog)

            if "qc" in report_json:
                ptr_qc = report_json["qc"]
                ptr["percent_reads_mapped"] = ptr_qc["percent_reads_mapped"]
                ptr["num_mapped_reads"] = ptr_qc["num_reads_mapped"]
                ptr["target_median_depth"] = ptr_qc["target_median_depth"]
                ptr["genome_median_depth"] = ptr_qc["genome_median_depth"]

            if "main_lineage" in report_json:

                tmp_ptr = {"lineage": [], "family": [], "rd": [], "fraction": [
                ], "spoligotype": [report_json.get("spoligotype", None)]}

                prv_lin = ""
                lineages = sorted(
                    report_json["lineage"], key=lambda x: x['lineage'])
                for lng in lineages[::-1]:
                    if not prv_lin.startswith(lng['lineage']):
                        tmp_ptr['fraction'].append(round(lng['fraction'], 2))
                        tmp_ptr['lineage'].append(lng['lineage'])
                        tmp_ptr['family'].append(lng['family'])
                        tmp_ptr['rd'].append(lng['rd'])
                    prv_lin = lng['lineage']

                is_mixed = len(tmp_ptr["lineage"]) > 1
                for k, v in tmp_ptr.items():
                    tmp_ptr[k] = '|'.join(map(str, v)).replace('None', '')
                tmp_ptr["is_mixed"] = is_mixed
                ptr['lineage'] = tmp_ptr

            if catalog.startswith('WHO'):
                fetch_variants_who(ptr, report_json, catalog,catch_resis_drg)
            elif catalog.startswith('TB-Profiler'):
                fetch_variants_tbdb(ptr, report_json, catalog)

            target_qcs = report_json["qc"]["target_qc"]
            deleted_genes = []

            for trgt in target_qcs:
                gene = trgt['target']
                if trgt['median_depth'] <= 5:
                    aff_drugs = gene_del_drugs(gene)
                    deleted_genes.append(gene)
                    if len(aff_drugs) > 0:
                        for drg, pred in aff_drugs:
                            ptr['variants'].append({
                                'catalog':catalog,
                                "drug": str(drg).upper(),
                                "gene_name": gene,
                                "change": "feature_ablation",
                                "effect": "feature_ablation",
                                "variant": f"{gene}_deletion",
                                "prediction": pred
                            })
                    else:
                        ptr["notes"].append(
                            f"{gene}: feature_ablation with unkown effect on resistance.")

            # unverified mutations
            ptr_unver = ptr['unverified_mutations']
            missing_locus = [convert_gene_locus(g) for g in deleted_genes]
            miss_ptr = report_json["qc"]["missing_positions"]
            for mis in miss_ptr:
                depth = mis['depth']
                pos = mis['pos']
                for annot in mis['annotation']:
                    locus = annot['gene']
                    # The gene is deleted no need to list related mutations as unverified
                    if locus in missing_locus:
                        continue

                    var = annot["variant"]
                    key = f"{gene}_{var}_{catalog}"
                    if key not in ptr_unver:
                        ptr_unver[key] = {
                            'catalog': catalog,
                            "depth": depth,
                            "pos": pos,
                            "drug": annot.get("drug", "").upper(),
                            "gene": convert_gene_locus(locus, True),
                            "locus": locus,
                            "variant": var,
                            "effect": annot.get("effect"),
                            "grade": annot.get("grade"),
                            "confidence": annot.get("confidence"),
                            "comment": annot.get("comment"),
                            "prediction": annot.get("prediction"),
                        }

    return report_dict


def fetch_variants_tbdb(json_ptr, report_json, catalog):
    sel_keys_base = ['pos', 'depth', 'freq', 'gene_id',
                     "gene_name", "change", "nucleotide_change", "protein_change"]
    sel_key_annot = ['drug', 'variant', 'effect', 'grade',
                     'prediction', 'comment_1', 'comment_2', 'footnote']
    conf_pred = {"Assoc w R": 'R', "Assoc w R - Interim": 'R',
                 "Uncertain significance": 'U', "Not assoc w R - Interim": 'S', "Not assoc w R": 'S'}
    
    ptr_var = json_ptr['variants']
    for rvr in report_json["dr_variants"]:
        base_info = {'catalog': catalog}
        for k in sel_keys_base:
            val = rvr.get(k, None)
            if k == 'freq' and val is not None:
                val = round(val, 2)
            
            base_info[k] = val

        for consq in rvr["consequences"]:
            for annot in consq["annotation"]:
                annot_content = {}
                for k in sel_key_annot:
                    val = annot.get(k, None)
                    if k == 'drug':
                        val = str(val).upper()
                    if k == 'effect':
                        val = rvr.get("type",None)
                    annot_content[k] = val

                annot_content['variant'] = annot['original_mutation']
                annot_content['comment_1'] = annot['comment']
                annot_content['grade'] = annot['confidence']
                annot_content['footnote'] = annot['source']

                if annot['type'] == "who_confidence":
                    annot_content['prediction'] = conf_pred.get(
                        annot['confidence'], None)
                elif annot['type'] == "drug_resistance":
                    annot_content['prediction'] = conf_pred.get(
                        annot['confidence'], 'R')
            
                if annot_content["prediction"] != 'S':
                    ptr_var.append({**base_info, **annot_content})

    for rvr in report_json["other_variants"]:

        if "annotation" not in rvr:
            continue

        base_info = {'catalog': catalog}
        for k in sel_keys_base:
            val = rvr.get(k, None)
            if k == 'freq' and val is not None:
                val = round(val, 2)
            base_info[k] = val

        for annot in rvr["annotation"]:

            annot_content = {}
            for k in sel_key_annot:
                val = annot.get(k, None)
                if k == 'drug':
                    val = str(val).upper()
                annot_content[k] = val

            annot_content['variant'] = annot['original_mutation']
            annot_content['comment_1'] = annot['comment']
            annot_content['grade'] = annot['confidence']
            annot_content['footnote'] = annot['source']

            if annot['type'] == "who_confidence":
                annot_content['prediction'] = conf_pred.get(
                    annot['confidence'], None)
            elif annot['type'] == "drug_resistance":
                annot_content['prediction'] = conf_pred.get(
                    annot['confidence'], 'R')
            if annot_content["prediction"] != 'S':
                ptr_var.append({**base_info, **annot_content})


def fetch_variants_who(json_ptr, report_json, catalog,save_resist_drug=False):
    sel_keys_base = ['pos', 'depth', 'freq', 'gene_id',
                     "gene_name", "change", "nucleotide_change", "protein_change"]
    sel_key_annot = ['drug', 'variant', 'effect', 'grade',
                     'prediction', 'comment_1', 'comment_2', 'footnote']

    ptr_var = json_ptr['variants']
    for rvr in report_json["dr_variants"]:
        base_info = {'catalog': catalog}
        for k in sel_keys_base:
            val = rvr.get(k, None)
            if k == 'freq' and val is not None:
                val = round(val, 2)
            base_info[k] = val

        for consq in rvr["consequences"]:
            for annot in consq["annotation"]:
                if annot["prediction"] in ['R', 'U']:
                    annot_content = {}
                    for k in sel_key_annot:
                        val = annot.get(k, None)
                        if k == 'drug':
                            val = str(val).upper()
                            if annot['prediction'] == 'R' and save_resist_drug and val not in json_ptr['resist_drug']:
                                json_ptr['resist_drug'].append(val)

                        annot_content[k] = val
                    ptr_var.append({**base_info, **annot_content})

    for rvr in report_json["other_variants"]:
        # ignore if there is no annotation:
        # TODO: would be good to investigate these mutations in the future.
        if "annotation" not in rvr:
            continue

        base_info = {'catalog': catalog}
        for k in sel_keys_base:
            val = rvr.get(k, None)
            if k == 'freq' and val is not None:
                val = round(val, 2)
            base_info[k] = val

        for annot in rvr["annotation"]:

            if annot['prediction'] != 'S':
                annot_content = {}
                for k in sel_key_annot:
                    val = annot.get(k, None)
                    if k == 'drug':
                        val = str(val).upper()
                    annot_content[k] = val
                ptr_var.append({**base_info, **annot_content})


def get_gene_del(region_qc):
    del_genes = {}
    for item in region_qc:
        gene = item['target']
        if item['median_depth'] <= 5:
            del_genes[gene] = gene_del_drugs(gene)
    return del_genes


def collect_cryptics(rep_dir, report_dict):
    cryp_out_lst = glob(os.path.join(rep_dir, "*cryptic.json"))

    for fj in cryp_out_lst:
        sid = os.path.basename(fj).split('.')[0]
        cat = os.path.basename(fj).split('.')[1]
        cat = cat.split('__')[-1]
        catalog = get_catalog_fullname('cryptic')

        if sid not in report_dict:
            report_dict[sid] = get_base_dict(sid)

        ptr = report_dict[sid]
        ptr['analyzed_catalogs'].append(catalog)
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
                        ptr_var.append({'catalog': catalog, "gene_name": mut['GENE'],
                                        'prediction': mut['PREDICTION'], 'drug': drg,
                                        "change": mut['MUTATION'],
                                        "footnote":cat})

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
    is_second_line_injectable_resistant = second_line_injectables.intersection(resistant_drugs)
    is_fluoroquinolone_resistant = fluoroquinolones.intersection(resistant_drugs)
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
    elif is_INH_resistant or is_RIF_resistant:
        return 'Pre-MDR'
    else:
        return 'other'


def generate(in_dir, prefix):
    report_dict = { "catalogs_repo": get_catalog_repo()}
    samples = {}
    collect_tbprofilers(in_dir, samples)
    collect_cryptics(in_dir, samples)

    report_dict['samples'] = list(samples.values())
    for i in range(len(report_dict['samples'])):
        sm = report_dict['samples'][i]
        sm['category_who'] = classify_tb_resistance(sm['resist_drug'])
        

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
