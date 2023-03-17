"""
This script gets the mean performances of multiple runs from combinations of 2 samples
"""

import os
import time
import pandas as pd
import sys


label_map = {
    "20230123_lb/": "lb",
}
append_map = {"20230123_lb/": []}
mageck_map = {"20230123_lb/": "20230123_lb/"}
model_folder = label_map.keys()

tb_rep_combs = ["rep1_rep2", "rep1_rep3_VPA", "rep2_rep3_VPA"]
fs_reps = [
    "rep5_rep8",
    "rep8_rep10",
    "rep5_rep11",
    "rep8_rep10",
    "rep8_rep11",
    "rep10_rep11",
]

prefix = [f"LDLvar_rep1_3_5_7_16.masked.{rc}" for rc in tb_rep_combs + fs_reps]
output_folder = "results/"
output_prefix = sys.argv[1]
for mf in model_folder:
    mmfs = append_map[mf]
    for p in prefix:
        os.system(
            f"python get_performance.py -g {output_folder}/{mf} {p} -m {mageck_map[mf]} &"
        )
        for mmf in mmfs:
            os.system(f"python get_performance.py -g {output_folder}/{mmf} {p} &")
time.sleep(20)
print("aggregating")
completed = False
trials = 0

while not completed and trials < 5:
    try:
        trials += 1
        dfs = []
        for mf in model_folder:
            print(f"Reading result of model {mf}")
            agg_metrics_dfs = []
            all_metrics_dfs = []
            for p in prefix:
                l = label_map[mf]
                df1 = pd.read_excel(
                    f"{output_folder}/{mf}/{p}.metrics.xlsx", index_col=0, header=[0, 1]
                )
                df1["run"] = p.split(".")[-1]
                for mmf in append_map[mf]:
                    df1_append = pd.read_excel(
                        f"{output_folder}/{mmf}/{p}.metrics.xlsx",
                        index_col=0,
                        header=[0, 1],
                    )
                    df1_append["run"] = p.split(".")[-1]
                    df1 = pd.concat([df1, df1_append])
                agg_metrics_dfs.append(df1)

                df2 = pd.read_excel(
                    f"{output_folder}/{mf}/{p}.all_metrics.xlsx",
                    index_col=0,
                    header=[0, 1],
                )
                df2["run"] = p.split(".")[-1]
                for mmf in append_map[mf]:
                    df2_append = pd.read_excel(
                        f"{output_folder}/{mmf}/{p}.all_metrics.xlsx",
                        index_col=0,
                        header=[0, 1],
                    )
                    df2_append["run"] = p.split(".")[-1]
                    df2 = pd.concat([df2, df2_append])
                all_metrics_dfs.append(df2)
        agg_metrics = pd.concat(agg_metrics_dfs)
        all_metrics = pd.concat(all_metrics_dfs)

        def write_mean_std(metrics: pd.DataFrame, writer_path: str):
            writer = pd.ExcelWriter(writer_path, engine="xlsxwriter")
            metrics.groupby(metrics.index).mean().to_excel(writer, sheet_name="mean")
            metrics.groupby(metrics.index).std().to_excel(writer, sheet_name="std")
            writer.close()
        agg_metrics.to_csv(f"{output_folder}/{output_prefix}_agg_metrics_comb2.csv")
        all_metrics.to_csv(f"{output_folder}/{output_prefix}_all_metrics_comb2.csv")
        completed = True

    except Exception as e:
        print(e)
        print("Reading file failed. Wait and try again.")
        time.sleep(10)
