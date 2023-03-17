import pandas as pd

samples = ["SubA_Q5", "SubB_Q5", "SubC_Q5", "SubC_Q5U", "SubD_Q5", "SubD_Q5U"]
reps = ["rep1","rep2","rep3", "WT_Ctrl"]
for s in samples:
    df_list = []
    for r in reps:
        tbl = pd.read_csv("CRISPResso_on_HepG2_LDLREndo_{}_{}_Endogenous_gDNA_R1/Quantification_window_nucleotide_percentage_table.txt".format(
            s, r), sep = "\t")
        tbl['replicate'] = r
        df_list.append(tbl)
    stats = pd.concat(df_list, axis=0)
    stats.to_csv('{}_qw_nt_perc.txt'.format(s), sep = "\t")

