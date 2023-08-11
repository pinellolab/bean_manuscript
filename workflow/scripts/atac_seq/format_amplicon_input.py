import pandas as pd

amplicon_info_tbl = pd.read_csv("resources/atac_seq/030123_ATACseq_info.csv")
guide_assignment_tbl = pd.read_excel(
    "resources/atac_seq/0223_LDLvar_indiv_summary.xlsx"
)
guide_assignment_tbl["guide_n"] = guide_assignment_tbl["individual gRNA name"].map(
    lambda s: s.split("g")[-1].split("_")[0]
)
guide_assignment_tbl["guide_name"] = (
    guide_assignment_tbl["Full target name"] + "_g" + guide_assignment_tbl["guide_n"]
)
guide_assignment_tbl["Full target name"] = guide_assignment_tbl[
    "Full target name"
].astype(str)
guide_assignment_tbl = guide_assignment_tbl.set_index("Full target name")
gRNA_seqs = pd.read_csv("resources/gRNA_info/LDLvar_gRNA_bean.csv")
amplicon_info_tbl = amplicon_info_tbl[
    ["gRNA", "Amplicon_fw (reference)", "Amplicon_rv (reference)"]
]
amplicon_info_tbl["gRNA"] = amplicon_info_tbl["gRNA"].astype(str)
amplicon_info_tbl["guide_name"] = guide_assignment_tbl.loc[
    amplicon_info_tbl.gRNA, "guide_name"
].reset_index(drop=True)


amplicon_info_tbl.columns = ["variant", "Amplicon_fw", "Amplicon_rv", "gRNA"]
# amplicon_info_tbl = pd.wide_to_long(
#     amplicon_info_tbl, stubnames="Amplicon_", i="variant", j="category", suffix="\\w+"
# ).reset_index()

amplicon_info_tbl = amplicon_info_tbl.merge(
    gRNA_seqs[["name", "sequence"]].rename(columns={"name": "gRNA"}),
    how="left",
    on="gRNA",
)
amplicon_info_tbl = amplicon_info_tbl[["variant", "Amplicon_fw", "sequence"]]
amplicon_info_tbl.columns = [
    "amplicon_name",
    "amplicon_seq",
    "guide_seq",
]
amplicon_info_tbl["amplicon_seq"] = amplicon_info_tbl["amplicon_seq"].map(
    lambda x: x.replace("N", " ").strip().replace(" ", "N")
)
amplicon_info_tbl.to_csv(
    "resources/atac_seq/030123_amplicon_info_CRISPRessoInput_fw_w_guide.tsv",
    sep="\t",
    index=False,
)

# tbl = amplicon_info_tbl[["AMPLICON_NAME", "AMPLICON_SEQUENCE"]].drop_duplicates()
# tbl.loc[tbl["AMPLICON_NAME"].map(lambda s: s.endswith("fw"))].to_csv(
#     "resources/atac_seq/030123_amplicon_info_CRISPRessoInput_fw_w_guide.tsv",
#     sep="\t",
#     index=False,
# )
