# -*- coding: utf-8 -*-
"""共享引用模块：PubMed 全量元数据引用串（终审版本）"""
import json
_FULL = json.load(open("code/citations_full.json"))

def cite_by_pmid(pmid):
    r = _FULL[pmid]
    au = r["authors"]
    title = r["title"].rstrip(".")
    end = "" if title.endswith("?") else "."
    pages = str(r.get("pages","")).rstrip("-.")
    doi = r.get("doi")
    doi_s = f" doi: {doi}." if doi else ""
    au = au.rstrip('.')
    return f"{au}. {title}{end} {r['journal']}. {r['year']};{r.get('vol') or r.get('volume','')}:{pages}.{doi_s} PMID: {pmid}."

# 已验证键 → PMID 映射（三稿共用）
KEY2PMID = {
    "randweaver2013": "24006913",
    "margiotta2014": "25338069",
    "pubchem2025": "39558165",
    "anastrozole_water": "33254902", "medaka_antiestrogen": "21963593",
    "ai_env_sludge": "42240727", "letrozole_zebrafish": "38452838", "ees_letrozole_tilapia": "33915448", "aromatase_fish_sensitivity": "31099392",
    "jumper2021": "34265844", "af2db2024": "37933859",
  # T1
  "sung2021": "33538338", "janjigian2021": "34102137", "gastric_ici2026": "42583108",
  "gastric_meta2025": "40299767", "hoadley2018": "29625048", "joly2003": "12942076",
  "ai_bone_faers2026": "42403084", "takeda2002": "12419242", "yadav2008": "20392224", "vanpuijenbroek2002": "11998548", "hu2021_oxc": "33446563", "vivian2017_toil": "28398314", "kommoss2017": "28159814", "lisowska2014": "24478986", "yoshihara2010": "20300634", "denkert2009": "19294737", "konstantinopoulos2010": "20547991", "mok2009": "19962670", "ferriss2012": "22348014", "zhang2022_scadv": "35196078", "keynote100": "31046082", "prima2019": "31562799", "tcgaov2011": "21720365", "simon2011": "27065756", "jung2019": "31537801", "gide2019": "30753825", "cho2020": "32879421", "solo2main": "28754483", "trogo2025": "39741180", "ctla4_2018": "29118008", "gulp2014": "24743597",
  "bindea2013": "24138885", "ssgsea2013": "23323831", "ccp2010": "20427518",
  "glmnet2010": "20808728", "timeroc2013": "24027076", "tide2018": "30127393",
  "tidepy2020": "32102694", "gdsc2013": "23180760", "onco2021": "34260682",
  "limma2015": "25605792", "geo2002": "11752295", "ooi2009": "19798449",
  "seurat2024": "37231261", "cellchat2021": "33597522", "stfk2022": "35510185",
  "l10002016": "28413689", "vina2021": "34278794", "pdb2002": "12037327",
  "nrp1_2026": "41311023", "sig2025": "39781724", "kim2018": "30013197",
  "barbie2009": "19847166",
  "ayers2017": "28650338",
  "tcga2014": "25079317",
  "cristescu2015": "25894828",
  "sohn2017": "28747339",
  "fuchs2018": "29543932",
  "kang2017": "28993052",
  "hamieh2019": "30918399",
  "miyake2021": "34069602",
  "bang2010": "20728210",
  "ies2025": "41081048",
  "sexton2020": "32894370",

  "goldman2020": "32444850", "iorio2016": "27397505",
  # T2
  "lheureux2019": "30910306", "perren2011": "22204725", "moore2026": "42705254",
  "gtex2020": "32913098",
  "dupont2011": "21654799",
  "zanconato2016": "27300434",
  "moore2018": "30345884",
  "bowtell2015": "26493647",
  "hamanishi2015": "26351349",
  "oxc2026": "42251166",
  "oxc_ccr2026": "40705052",

  "ovtme2025": "40369674", "hamidi2018": "30002479", "panciera2017": "28951564",
  "piezo2026": "42836163", "lamina2026": "42467792", "bonome2008": "18593951",
  "scissor2022": "34764492", "ces1_2026": "42870294", "ov_ici2026": "42865715",
  "olaparib2026": "41861615",
  # T3
  "burstein2021": "34324367", "burstein2022": "35536015", "rachner2011": "21450337",
  "osteopath2026": "42833929", "romosozumab2026": "42274886", "ai_faers2026": "42541454",
  "cyp19_2027": "42735742", "hopkins2008": "18936753", "nogales2022": "34895945",
  "bate2009": "19358225", "openfda2016": "26644398", "chembl2024": "37933841",
  "swiss2019": "31106366", "genecards2016": "27322403", "string2023": "36370105",
  "wgcna2008": "19114008", "cp2021": "34557778", "boneniche2026": "42729472",
  "vps35_2021": "33997170", "retromer2023": "37350516", "sfrp4_2022": "35182301",
  "m6a_2024": "38658320",
  "boyce2007": "17634140",
  "ai_bone2025": "41568387",
  "coleman2020": "32801018",
  "hadji2025": "40726588",
  "seaman2012": "23148298",
  "baron2013": "23389618",

}

def refs_for(keys, extra=None):
    out = {}
    for k in keys:
        if k in KEY2PMID:
            out[k] = cite_by_pmid(KEY2PMID[k])
    if extra:
        out.update(extra)
    return out
