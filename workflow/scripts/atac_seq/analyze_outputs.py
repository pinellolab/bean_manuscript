import os
import pandas as pd

reps = [1, 2, 3, "C"]
subs = ["Ser", "SS"]
libraries = ["ATAC_Seq", "GDNA"]
results = []
allele_info = pd.read_excel(
    "resources/LDLvar/20221013_LDLvar_simpleZscores_credset.xlsx", sheet_name=1
)[["labels", "A1", "A2", "HepG2 alternate allele fraction"]]
for rep in reps:
    for sub in subs:
        locus = f"{rep}-{sub}_20-locus"
        os.system(f"python scripts/atac_seq/analyze_output.py {locus}")
        result = pd.read_csv(f"results/atac_seq/{locus}_edits.csv")
        results.append(result)
res_tbl = pd.concat(results)
res_tbl.merge(allele_info, left_on="amplicon", right_on="labels").to_csv(
    "results/atac_seq/all_samples_edits.csv"
)
