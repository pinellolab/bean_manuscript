import os
import sys
import multiprocessing
import numpy as np
import pandas as pd
from tqdm import tqdm
import beret as be
import swifter

from multiprocessing import set_start_method

allele_freq_tbl_prefix = "Alleles_frequency_table_around_"
prefix = sys.argv[1]
def convert_crispresso_to_Allele_row(row):
    allele = be.mp._supporting_fn._get_allele_from_alignment(row.Reference_Sequence,
        row.Aligned_Sequence,
        offset = 0,
        strand = 1,
        start_pos = 0,
        end_pos = np.inf,
        )
    return(allele)

def convert_crispresso_to_Allele(crispresso_allele_tbl: pd.DataFrame, guide, count_column_name = "sample"):
    alleles = crispresso_allele_tbl.apply(convert_crispresso_to_Allele_row, axis=1)
    allele_count_tbl = pd.DataFrame({"guide":guide, "allele":alleles, count_column_name:crispresso_allele_tbl["#Reads"]})
    allele_count_tbl = allele_count_tbl.groupby(['guide', 'allele']).sum().reset_index()
    return allele_count_tbl

def convert_sample_to_allele_df(sample_name = "sample"):
    crispresso_dir = "{}/CRISPResso_on_HepG2_LDLREndo_{}_Endogenous_gDNA_R1".format(prefix, sample_name)
    allele_tbl_paths = [filename for filename in os.listdir(crispresso_dir) if filename.startswith(allele_freq_tbl_prefix)]
    guides = [filename.split(allele_freq_tbl_prefix)[-1].split(".txt")[0] for filename in allele_tbl_paths]
    count_tbls = []
    for allele_tbl_path, guide in zip(allele_tbl_paths, guides):
        count_tbls.append(convert_crispresso_to_Allele(pd.read_table("{}/{}".format(crispresso_dir, allele_tbl_path)),
            guide,
            sample_name))
    allele_count_tbl = pd.concat(count_tbls, axis=0).reset_index(drop=True)
    return(allele_count_tbl)

if __name__ == '__main__':
    set_start_method("spawn", force=True)
    processes = []
    sublibrary_experiments = {"A":["Q5"],
        "B":["Q5"],
        "C":["Q5", "Q5U"],
        "D":["Q5", "Q5U"]}
    replicates = ["rep1", "rep2", "rep3", "WT_Ctrl"]

    samples = []
    for sublibrary in ['A', 'B', 'C', 'D']:
        for exp in sublibrary_experiments[sublibrary]:
            for rep in replicates:
                samples.append("Sub{}_{}_{}".format(sublibrary, exp, rep))

    with multiprocessing.Pool(30) as p:
        allele_dfs = p.map(convert_sample_to_allele_df, samples)
    
    
    # for i, sample in enumerate(samples):
    #     allele_df = allele_dfs[i]
    #     guide_counts = allele_df.groupby('guide')[sample].sum()
    #     guides = guide_counts.index.tolist()
    #     screen = be.ReporterScreen(X = np.reshape(guide_counts.to_numpy(), (-1,1)),
    #         uns={"allele_counts":allele_df},
    #         guides=pd.DataFrame({"name":guides}),
    #         condit=pd.DataFrame(index=[sample]))
    #     screen.write("{}/Sub{}.h5ad".format(prefix, sample))

    i=0
    for sublibrary in tqdm(['A', 'B', 'C', 'D'], desc = "sublibrary"):
        
        for exp in sublibrary_experiments[sublibrary]:
            screens = []
            for rep in tqdm(replicates, desc = "replicates", leave = False):
                allele_df = allele_dfs[i]
                sample = samples[i]
                guide_counts = allele_df.groupby('guide')[sample].sum()
                guides = guide_counts.index.tolist()
                screen = be.ReporterScreen(X = np.reshape(guide_counts.to_numpy(), (-1,1)),
                    uns={"allele_counts":allele_df},
                    guides=pd.DataFrame({"name":guides}),
                    condit=pd.DataFrame(index=[sample]))
                screens.append(screen)
                i+=1
            sample_screen = be.concat(screens)
            sample_screen.write("{}/Sub{}_{}.h5ad".format(prefix, sublibrary, exp))

