import pandas as pd

samples = ["SubA_Q5", "SubB_Q5", "SubC_Q5", "SubC_Q5U", "SubD_Q5", "SubD_Q5U"]
reps = ["rep1","rep2","rep3", "WT_Ctrl"]
df_list = []
label_list = []
for s in samples:
    for r in reps:
        tbl = pd.read_csv("CRISPResso_on_HepG2_LDLREndo_{}_{}_Endogenous_gDNA_R1/CRISPResso_mapping_statistics.txt".format(
            s, r), sep = "\t")
        df_list.append(tbl.T)
        label_list.append("{}_{}".format(r, s))
stats = pd.concat(df_list, axis=1)
stats.columns=label_list
stats.to_csv('mapping_stats.txt', sep = "\t")

