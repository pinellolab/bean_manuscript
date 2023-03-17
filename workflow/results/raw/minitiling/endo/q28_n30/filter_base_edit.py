import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

base_to_rowidx = {"A":0, "C":1, "G":2, "T":3}
def test_base_greater(treat_folder, ctrl_folder):
    quant_tbl_file = "Quantification_window_nucleotide_frequency_table.txt"
    tbl =pd.read_table("{}/{}".format(treat_folder, quant_tbl_file), index_col = 0).iloc[:4,:]
    refseq = tbl.columns.map(lambda s: s.split(".")[0]).tolist()
    count1 = pd.read_table("{}/{}".format(treat_folder, quant_tbl_file), index_col = 0).values
    count2 = pd.read_table("{}/{}".format(ctrl_folder, quant_tbl_file), index_col = 0).values
    
    odds_ratio_tbl = pd.DataFrame().reindex_like(tbl)
    p_value_tbl = pd.DataFrame().reindex_like(tbl)
    for i in range(len(tbl.T)):
        ref_row = base_to_rowidx[refseq[i]]
        for j in range(4):
            if j == ref_row: 
                odds_ratio_tbl.iloc[j, i] = np.nan
                p_value_tbl.iloc[j, i] = np.nan     
            else:
                fisher_tbl = [[count1[j, i], count1[:, i].sum() - count1[j, i]],
                [count2[j, i], count2[:, i].sum() - count2[j, i]]]
                odds_ratio_tbl.iloc[j, i], p_value_tbl.iloc[j, i] = fisher_exact(
                    fisher_tbl, alternative = "greater")
    q_bonf_tbl = p_value_tbl * np.count_nonzero(~np.isnan(p_value_tbl.values))
    q_bonf_tbl[q_bonf_tbl > 1] = 1
    return q_bonf_tbl, odds_ratio_tbl

def _get_mask_base_substitution(treat_folder, ctrl_folder, q_thres = 0.001, OR_thres = 5):
    q, odds_ratio = test_base_greater(treat_folder, ctrl_folder)
    return(((q < q_thres) & (odds_ratio > OR_thres)), q, odds_ratio)

def filter_base_substitution(treat_folder, ctrl_folder, **kwargs):
    quant_perc_tbl = "Quantification_window_nucleotide_percentage_table.txt"
    perc_tbl = pd.read_table("{}/{}".format(treat_folder, quant_perc_tbl), index_col = 0).iloc[:4,:] * 100
    
    mask, q, odds_ratio = _get_mask_base_substitution(treat_folder, ctrl_folder, **kwargs)
    perc_tbl[~mask] = np.nan
    return(perc_tbl, q, odds_ratio)