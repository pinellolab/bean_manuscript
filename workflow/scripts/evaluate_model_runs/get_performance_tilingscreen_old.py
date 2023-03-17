"""Get performance of tiling screen results."""
import argparse
import re
import sys

import anbe_tool.evaluate as ae
import anbe_tool.translate as at
import bean as be
import numpy as np
import pandas as pd


def decode_control_edits(edit_str: str, guide_len: int) -> str:
    """For edit with unique IDs for negative control (ex. CBE_CONTROL_100_pos!13:A>G),
    format the edit string
    """
    # TODO: match position by transforming from rel_pos based on reporter to guide
    assert re.fullmatch(r"[\w_]*!\d+:[A-Z*-]>[A-Z*-]", edit_str), edit_str
    uid, edit = edit_str.split("!")
    pos, transition = edit.split(":")
    guide_pos = int(pos) - (32 - 6 - guide_len)
    return f"{uid}_{guide_pos}:{transition}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get performance metric and plot.")
    parser.add_argument(
        "model_folder", type=str, help="relative path to model run folder."
    )
    parser.add_argument(
        "prefix",
        type=str,
        help="file description for the run. (ex. LDLvar_rep1_5_7_16.masked)",
    )
    parser.add_argument(
        "-m",
        "--mageck_model_folder",
        type=str,
        default=None,
        help="relative path to Mageck run folder.",
    )
    parser.add_argument("--plot", "-g", action="store_true", help="plot the table")
    parser.add_argument(
        "--print_res", "-p", action="store_true", help="Print out the table"
    )

    args = parser.parse_args()
    model_folder = args.model_folder
    prefix = args.prefix
    plot = args.plot
    mageck_fdr_col = "wald-fdr"

    model_ids = ["C16", "C18", "C19"]
    res_tbls = []
    for model_id in model_ids:
        result_path = (
            f"{model_folder}/{prefix}.model{model_id}.CRISPRbean_element_result.csv"
        )
        try:
            tbl = pd.read_csv(result_path)
            if len(res_tbls) == 0:
                res_tbls.append(tbl.iloc[:, :10])
            res_tbls.append(tbl.iloc[:, 10:].add_suffix(f"_{model_id}"))
        except FileNotFoundError:
            print(f"couldn't find {result_path}: ignoring")
            pass

    if len(res_tbls) == 0:
        print(f"No result found for {model_folder}/{prefix}")
        exit(1)

    fit_df = pd.concat(res_tbls, axis=1)
    guide_info_df = pd.read_csv(
        "~/projects/ANBE/data/gRNA_info/LDLRCDS_gRNA_beret.csv"
    ).set_index(
        "name"
    )  # TODO allow users to feed this in?
    guide_lens = guide_info_df.loc[
        guide_info_df.Region.isin(["CBE control", "ABE control"]), "guide_len"
    ]
    fit_df["edit"] = fit_df["edit"].map(
        lambda s: decode_control_edits(s, guide_lens[s.split("!")[0]])
        if "CONTROL" in s
        else s
    )
    fit_df.to_csv(f"{model_folder}/fit_df_{prefix}.csv")

    labels = [f"mu_z_adj_{model_id}" for model_id in model_ids]

    for i, mageck_outcome_type in enumerate(["allEdited", "behive"]):
        mageck_modes = ["topbot", "topbot_var", "sort", "sort_var"]
        mageck_labels = []
        mageck_results = None
        # for mageck_outcome_type in ["behive"]:
        if not args.mageck_model_folder is None:
            # Reading MAGeCK MLE
            mageck_prefix = f"/PHShome/jr1025/projects/ANBE/mageck_results/{args.mageck_model_folder}/{prefix}.target_{mageck_outcome_type}/"
            for mageck_mode in mageck_modes:
                tbl = pd.read_csv(
                    f"{mageck_prefix}/{mageck_mode}.gene_summary.txt",
                    sep="\t",
                    index_col=0,
                )
                tbl.columns = tbl.columns.map(lambda s: s.split("|")[-1])
                tbl = tbl.add_suffix(f"_{mageck_outcome_type}.{mageck_mode}")
                if mageck_outcome_type == "behive":
                    pass
                    # # TODO: remove this after re-running mageck with correct edit assignment
                    # tbl.index = tbl.index.map(
                    #     lambda s: f"A{s}"
                    #     if re.fullmatch(r"\d+:[A-Z-*]>[A-Z-*]", s)
                    #     and int(s.split(":")[0].split("!")[-1]) < 1e4
                    #     else s
                    # )
                else:
                    # TODO: remove this after re-running mageck with correct edit assignment
                    def fix_control_pos(s):
                        if "CONTROL" in s:
                            uid_pos, transition = s.split(":")
                            uid, pos = uid_pos.rsplit("_", 1)
                            return f"{uid}_{int(pos)+1}:{transition}"
                        else:
                            return s
                    tbl.index = tbl.index.map(fix_control_pos)
                if mageck_results is None:
                    mageck_results = tbl
                else:
                    mageck_results = mageck_results.join(tbl, how="outer")

                mageck_labels.append(f"{mageck_outcome_type}.{mageck_mode}")
            result_combined_common = fit_df.merge(
                mageck_results, left_on="edit", right_on="Gene", how="inner"
            )

        else:
            result_combined_common = fit_df
        if i == 0:
            result_combined_all = fit_df.merge(
                mageck_results.reset_index(),
                left_on="edit",
                right_on="Gene",
                how="outer",
            )
        else:
            result_combined_all = result_combined_all.merge(
                mageck_results.reset_index(),
                left_on="edit",
                right_on="Gene",
                how="outer",
            )
        result_combined_all.to_csv(f"{model_folder}/results_merged_{prefix}.csv")
        result_df = result_combined_common

        list_zs = [
            result_df[f"mu_z_adj_{model_id}"].values for model_id in model_ids
        ] + [result_df[f"z_{mageck_label}"].values for mageck_label in mageck_labels]
        list_fdrs_inc = [
            result_df[f"fdr_inc_adj_{model_id}"].values for model_id in model_ids
        ] + [result_df[f"z_{mageck_label}"].values for mageck_label in mageck_labels]
        list_fdrs_dec = [
            result_df[f"fdr_dec_adj_{model_id}"].values for model_id in model_ids
        ] + [
            result_df[f"{mageck_fdr_col}_{mageck_label}"].values
            for mageck_label in mageck_labels
        ]
        splicing_idx = np.where(result_df.group == "splicing")[0]
        negctrl_idx = np.where(result_df.group == "negctrl")[0]
        syn_idx = np.where(result_df.group == "syn")[0]
        metric_column_labels = [
            f"Negative Control (n={len(negctrl_idx)})",
            f"Synonymous (n={len(syn_idx)})",
        ]
        labels_use = labels + mageck_labels
        precision = pd.concat(
            [
                ae.get_metrics(
                    ae._get_precision_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    negctrl_idx,
                ).round(3),
                ae.get_metrics(
                    ae._get_precision_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    syn_idx,
                ).round(3),
            ],
            axis=1,
        )
        precision.columns = ["negctrl", "syn"]

        recall = pd.concat(
            [
                ae.get_metrics(
                    ae._get_recall_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    negctrl_idx,
                ).round(3),
                ae.get_metrics(
                    ae._get_recall_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    syn_idx,
                ).round(3),
            ],
            axis=1,
        )
        recall.columns = ["negctrl", "syn"]

        fdr_f1 = pd.concat(
            [
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    negctrl_idx,
                ).round(3),
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    syn_idx,
                ).round(3),
            ],
            axis=1,
        )
        fdr_f1.columns = ["negctrl", "syn"]

        fdr_f01 = pd.concat(
            [
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    negctrl_idx,
                    beta=0.1,
                ).round(3),
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    syn_idx,
                    beta=0.1,
                ).round(3),
            ],
            axis=1,
        )
        fdr_f01.columns = ["negctrl", "syn"]
        aurocs = pd.concat(
            [
                ae.get_metrics(
                    ae._get_auroc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    negctrl_idx,
                    beta=0.1,
                ).round(3),
                ae.get_metrics(
                    ae._get_auroc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    syn_idx,
                    beta=0.1,
                ).round(3),
            ],
            axis=1,
        )
        aurocs.columns = ["negctrl", "syn"]

        auprcs = pd.concat(
            [
                ae.get_metrics(
                    ae._get_auprc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    negctrl_idx,
                    beta=0.1,
                ).round(3),
                ae.get_metrics(
                    ae._get_auprc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    splicing_idx,
                    syn_idx,
                    beta=0.1,
                ).round(3),
            ],
            axis=1,
        )
        auprcs.columns = ["negctrl", "syn"]

        perf_dict = {
            "Precision": precision,
            "Recall": recall,
            "F1": fdr_f1,
            "F0.1": fdr_f01,
            "AUROC": aurocs,
            "AUPRC": auprcs,
        }
        perf_df = pd.concat(perf_dict.values(), axis=1, keys=perf_dict.keys())
        if args.print_res:
            print(perf_df.to_csv(sep="\t"))

        for label, df in zip(
            ["all", "Precision", "Recall", "F1", "F0.1", "auroc", "auprc"],
            [perf_df, precision, recall, fdr_f1, fdr_f01, aurocs, auprcs],
        ):
            df.to_csv(f"{model_folder}/{prefix}.{label}.tsv", sep="\t")
        if plot:
            print(perf_df)
            # avg_d = {
            #     "Precision": get_average_metric(precision),
            #     "Recall": get_average_metric(recall),
            #     "F1": get_average_metric(fdr_f1),
            #     "F0.1": get_average_metric(fdr_f01),
            #     "AUROC": get_average_metric(aurocs),
            #     "AUPRC": get_average_metric(auprcs),
            # }
            # avg_df = pd.concat(avg_d.values(), axis=1, keys=avg_d.keys())
            # avg_df.style.set_properties(**{"width": "px"}).background_gradient().to_excel(
            #     f"{model_folder}/{prefix}.metrics.xlsx"
            # )
            outfile_path = f"{model_folder}/{prefix}.splicing.{mageck_outcome_type}.all_metrics.xlsx"
            perf_df.style.set_properties(
                **{"width": "px"}
            ).background_gradient().to_excel(outfile_path)
            print(f"Result saved into {outfile_path}")

        # Metric for pathogenic/benign variants
        result_df_coding = at.get_dms_df(result_df.loc[result_df.coding == "coding"])
        if i == 0:
            result_combined_all = at.get_dms_df(result_combined_all, join_how="left")
        else:
            result_combined_all.to_csv(f"{model_folder}/{prefix}.all_scores.csv")
            print(f"Scores saved at {model_folder}/{prefix}.all_scores.csv")

        list_zs = [
            result_df_coding[f"mu_z_adj_{model_id}"].values for model_id in model_ids
        ] + [
            result_df_coding[f"z_{mageck_label}"].values
            for mageck_label in mageck_labels
        ]
        list_fdrs_inc = [
            result_df_coding[f"fdr_inc_adj_{model_id}"].values for model_id in model_ids
        ] + [
            result_df_coding[f"z_{mageck_label}"].values
            for mageck_label in mageck_labels
        ]
        list_fdrs_dec = [
            result_df_coding[f"fdr_dec_adj_{model_id}"].values for model_id in model_ids
        ] + [
            result_df_coding[f"{mageck_fdr_col}_{mageck_label}"].values
            for mageck_label in mageck_labels
        ]
        lpidx = np.where(
            result_df_coding.clinvar_annot_2 == "Pathogenic/Likely_Pathogenic"
        )[0]
        pidx = np.where(result_df_coding.clinvar_annot_3 == "Pathogenic")[0]
        bidx = np.where(result_df_coding.clinvar_annot_2 == "Benign/Likely_Benign")[0]
        metric_column_labels = [
            f"Likely Pathogenic/Pathogenic (n={len(lpidx)})",
            f"Pathogenic (n={len(pidx)})",
        ]
        precision = pd.concat(
            [
                ae.get_metrics(
                    ae._get_precision_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                ).round(3),
                ae.get_metrics(
                    ae._get_precision_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                ).round(3),
            ],
            axis=1,
        )
        precision.columns = metric_column_labels

        recall = pd.concat(
            [
                ae.get_metrics(
                    ae._get_recall_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                ).round(3),
                ae.get_metrics(
                    ae._get_recall_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                ).round(3),
            ],
            axis=1,
        )
        recall.columns = metric_column_labels

        fdr_f1 = pd.concat(
            [
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                ).round(3),
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                ).round(3),
            ],
            axis=1,
        )
        fdr_f1.columns = metric_column_labels

        fdr_f01 = pd.concat(
            [
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                    beta=0.1,
                ).round(3),
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                    beta=0.1,
                ).round(3),
            ],
            axis=1,
        )
        fdr_f01.columns = metric_column_labels
        aurocs = pd.concat(
            [
                ae.get_metrics(
                    ae._get_auroc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                    beta=0.1,
                ).round(3),
                ae.get_metrics(
                    ae._get_auroc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                    beta=0.1,
                ).round(3),
            ],
            axis=1,
        )
        aurocs.columns = metric_column_labels

        auprcs = pd.concat(
            [
                ae.get_metrics(
                    ae._get_auprc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                    beta=0.1,
                ).round(3),
                ae.get_metrics(
                    ae._get_auprc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                    beta=0.1,
                ).round(3),
            ],
            axis=1,
        )
        auprcs.columns = metric_column_labels

        perf_dict = {
            "Precision": precision,
            "Recall": recall,
            "F1": fdr_f1,
            "F0.1": fdr_f01,
            "AUROC": aurocs,
            "AUPRC": auprcs,
        }
        perf_df = pd.concat(perf_dict.values(), axis=1, keys=perf_dict.keys())
        if args.print_res:
            print(perf_df.to_csv(sep="\t"))

        for label, df in zip(
            ["all", "Precision", "Recall", "F1", "F0.1", "auroc", "auprc"],
            [perf_df, precision, recall, fdr_f1, fdr_f01, aurocs, auprcs],
        ):
            df.to_csv(f"{model_folder}/{prefix}.{label}.tsv", sep="\t")
        if plot:
            print(perf_df)
            # avg_d = {
            #     "Precision": get_average_metric(precision),
            #     "Recall": get_average_metric(recall),
            #     "F1": get_average_metric(fdr_f1),
            #     "F0.1": get_average_metric(fdr_f01),
            #     "AUROC": get_average_metric(aurocs),
            #     "AUPRC": get_average_metric(auprcs),
            # }
            # avg_df = pd.concat(avg_d.values(), axis=1, keys=avg_d.keys())
            # avg_df.style.set_properties(**{"width": "px"}).background_gradient().to_excel(
            #     f"{model_folder}/{prefix}.metrics.xlsx"
            # )
            outfile_path = f"{model_folder}/{prefix}.clinvar.{mageck_outcome_type}.all_metrics.xlsx"
            perf_df.style.set_properties(
                **{"width": "px"}
            ).background_gradient().to_excel(outfile_path)
            print(f"Result saved into {outfile_path}")

        # Metric for pathogenic/benign variants
        result_df_clinvar = at.get_dms_df(
            result_df.loc[result_df.coding == "coding"], join_how="outer"
        )
        list_zs = [
            result_df_clinvar[f"mu_z_adj_{model_id}"].values for model_id in model_ids
        ] + [
            result_df_clinvar[f"z_{mageck_label}"].values
            for mageck_label in mageck_labels
        ]
        list_fdrs_inc = [
            result_df_clinvar[f"fdr_inc_adj_{model_id}"].values
            for model_id in model_ids
        ] + [
            result_df_clinvar[f"z_{mageck_label}"].values
            for mageck_label in mageck_labels
        ]
        list_fdrs_dec = [
            result_df_clinvar[f"fdr_dec_adj_{model_id}"].values
            for model_id in model_ids
        ] + [
            result_df_clinvar[f"{mageck_fdr_col}_{mageck_label}"].values
            for mageck_label in mageck_labels
        ]
        lpidx = np.where(
            result_df_clinvar.clinvar_annot_2 == "Pathogenic/Likely_Pathogenic"
        )[0]
        pidx = np.where(result_df_clinvar.clinvar_annot_3 == "Pathogenic")[0]
        bidx = np.where(result_df_clinvar.clinvar_annot_2 == "Benign/Likely_Benign")[0]
        metric_column_labels = [
            f"Likely Pathogenic/Pathogenic (n={len(lpidx)})",
            f"Pathogenic (n={len(pidx)})",
        ]
        precision = pd.concat(
            [
                ae.get_metrics(
                    ae._get_precision_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                ).round(3),
                ae.get_metrics(
                    ae._get_precision_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                ).round(3),
            ],
            axis=1,
        )
        precision.columns = metric_column_labels

        recall = pd.concat(
            [
                ae.get_metrics(
                    ae._get_recall_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                ).round(3),
                ae.get_metrics(
                    ae._get_recall_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                ).round(3),
            ],
            axis=1,
        )
        recall.columns = metric_column_labels

        fdr_f1 = pd.concat(
            [
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                ).round(3),
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                ).round(3),
            ],
            axis=1,
        )
        fdr_f1.columns = metric_column_labels

        fdr_f01 = pd.concat(
            [
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                    beta=0.1,
                ).round(3),
                ae.get_metrics(
                    ae._get_fdr_f_score_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                    beta=0.1,
                ).round(3),
            ],
            axis=1,
        )
        fdr_f01.columns = metric_column_labels
        aurocs = pd.concat(
            [
                ae.get_metrics(
                    ae._get_auroc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                    beta=0.1,
                ).round(3),
                ae.get_metrics(
                    ae._get_auroc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                    beta=0.1,
                ).round(3),
            ],
            axis=1,
        )
        aurocs.columns = metric_column_labels

        auprcs = pd.concat(
            [
                ae.get_metrics(
                    ae._get_auprc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    lpidx,
                    bidx,
                    beta=0.1,
                ).round(3),
                ae.get_metrics(
                    ae._get_auprc_dec,
                    list_zs,
                    list_fdrs_dec,
                    labels_use,
                    pidx,
                    bidx,
                    beta=0.1,
                ).round(3),
            ],
            axis=1,
        )
        auprcs.columns = metric_column_labels

        perf_dict = {
            "Precision": precision,
            "Recall": recall,
            "F1": fdr_f1,
            "F0.1": fdr_f01,
            "AUROC": aurocs,
            "AUPRC": auprcs,
        }
        perf_df = pd.concat(perf_dict.values(), axis=1, keys=perf_dict.keys())
        if args.print_res:
            print(perf_df.to_csv(sep="\t"))

        for label, df in zip(
            ["all", "Precision", "Recall", "F1", "F0.1", "auroc", "auprc"],
            [perf_df, precision, recall, fdr_f1, fdr_f01, aurocs, auprcs],
        ):
            df.to_csv(f"{model_folder}/{prefix}.{label}.tsv", sep="\t")
        outfile_path = f"{model_folder}/{prefix}.clinvar_all.{mageck_outcome_type}.all_metrics.xlsx"
        perf_df.style.set_properties(**{"width": "px"}).background_gradient().to_excel(
            outfile_path
        )
        print(f"Result saved into {outfile_path}")

        if plot:
            pass
            # Draw PRC curve

            # Draw Precision & Recall
