import sys
import pandas as pd
from Bio.Seq import Seq

"""
1. Get substitution output
2. Discard ref/alt allele and report it
"""
locus_name = sys.argv[1]  # "1-SS_20-locus"
window_start_pos = 2
window_end_pos = 7
gRNA_tbl_path = "resources/gRNA_info/LDLvar_gRNA_bean.csv"
amplicon_path = "resources/atac_seq/030123_amplicon_info_CRISPRessoInput_fw.tsv"
guide_assignment_path = "resources/atac_seq/0223_LDLvar_indiv_summary.xlsx"
amplicon_info_path = "_"

gdna_crispresso_prefix = f"results/atac_seq/CRISPRessoPooled_on_{locus_name}_GDNA_R1/"
atac_crispresso_prefix = (
    f"results/atac_seq/CRISPRessoPooled_on_{locus_name}_ATAC_Seq_R1/"
)


amplicon_info_tbl = pd.read_csv(amplicon_info_path)
sample_prefixes_dict = {
    "gdna": gdna_crispresso_prefix,
    "atac": atac_crispresso_prefix,
}

gRNA_tbl = pd.read_csv(gRNA_tbl_path, index_col=0).set_index("name")
guide_assignment = pd.read_excel(guide_assignment_path)
guide_assignment_tbl = guide_assignment[
    ["Full target name", "individual gRNA name"]
].copy()
guide_assignment_tbl["guide_n"] = guide_assignment_tbl["individual gRNA name"].map(
    lambda s: s.split("g")[-1].split("_")[0]
)
guide_assignment_tbl["guide_name"] = (
    guide_assignment_tbl["Full target name"] + "_g" + guide_assignment_tbl["guide_n"]
)
guide_assignment_tbl = guide_assignment_tbl.set_index("Full target name")

amplicon_tbl = pd.read_csv(amplicon_path, sep="\t")
amplicon_tbl["gRNA"] = amplicon_tbl["AMPLICON_NAME"].map(lambda s: s.rsplit("_", 1)[0])
amplicon_tbl = (
    amplicon_tbl[["gRNA", "AMPLICON_SEQUENCE"]]
    .merge(
        amplicon_info_tbl[["gRNA", "Varpos_fw", "Varnt_fw (reference)"]],
        on="gRNA",
        how="left",
    )
    .set_index("gRNA")
)
amplicon_tbl.columns = ["sequence", "pos", "ref_base"]
amplicons = amplicon_tbl.index.tolist()


def get_editing_window_allele(var_pos, var_refbase, nt_freq_tbl):
    edits_df = nt_freq_tbl.iloc[:, [var_pos]].T
    edits_df.columns = ["ref"] + nt_freq_tbl.index.tolist()[1:]
    return edits_df


no_files_output = open("results/atac_seq/no_crispresso_output.csv", "w")
edits_df_amplicons = []
for amplicon in amplicons:
    offset = 4 if amplicon == "rs771555783_Maj_ABE_309" else 0
    edits_df_dict = {}
    for sample_id, sample_prefix in sample_prefixes_dict.items():
        nt_freq_tbl_path = f"{sample_prefix}/CRISPResso_on_{amplicon}_fw/results_atac_seq_.Nucleotide_frequency_table.txt"
        try:
            nt_freq_tbl = pd.read_csv(
                nt_freq_tbl_path, sep="\t", header=None, index_col=0
            )
        except FileNotFoundError as exc:
            print(
                f"WARNING: {nt_freq_tbl_path} not found in {locus_name}. Ignoring the file..."
            )
            no_files_output.write(f"{nt_freq_tbl_path},{locus_name}")

        edits_df = get_editing_window_allele(
            amplicon_tbl.loc[amplicon, "pos"] - 4 + offset,
            amplicon_tbl.loc[amplicon, "ref_base"],
            nt_freq_tbl,
        )
        edits_df["locus"] = locus_name
        edits_df["sample"] = sample_id
        edits_df_dict[sample_id] = edits_df
    edits_df_cat = pd.concat(edits_df_dict.values())
    edits_df_cat["amplicon"] = amplicon
    edits_df_amplicons.append(edits_df_cat)
pd.concat(edits_df_amplicons).to_csv(f"results/atac_seq/{locus_name}_edits.csv")
no_files_output.close()
