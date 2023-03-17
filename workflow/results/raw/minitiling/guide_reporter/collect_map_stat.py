import pandas as pd

samples = ["SubA", "SubB", "SubC", "SubD"]
reps = [1,2,3]
lib_to_enz = {"SubA":["Q5"], "SubB":["Q5"], "SubC":["Q5", "Q5U"], "SubD":["Q5", "Q5U"]}
df_list = []
label_list = []
for lib in samples:
    for r in reps:
        for enz in lib_to_enz[lib]:
            tbl = pd.read_csv("{lib}_rawseq/rep{r}_{lib}_{enz}/crisprep_count_rep{r}_{lib}_{enz}/mapping_stats.txt".format(
                lib=lib, r=r, enz=enz), sep = "\t", comment="R", header=None, index_col = 0)
            df_list.append(tbl.T)
            label_list.append("rep{}_{}_{}".format(r, lib, enz))
stats = pd.concat(df_list)
stats.index=label_list
stats.to_csv('mapping_stats_rawseq.txt', sep = "\t")

