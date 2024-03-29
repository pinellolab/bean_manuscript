import argparse
import bean as be
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        prog="allele_filter",
        description="Filter alleles based on edit position in spacer and frequency across samples.",
    )
    parser.add_argument(
        "h5ad_path",
        help="Input ReporterScreen file of which allele will be filtered out.",
    )
    parser.add_argument(
        "output-prefix",
        help="Output prefix for log and ReporterScreen file with allele assignment",
    )
    parser.add_argument(
        "--plasmid_path",
        "-p",
        help="Plasmid ReporterScreen object path. If provided, alleles are filtered based on if a nucleotide edit is more significantly enriched in sample compared to the plasmid data. Negative control data where no edit is expected can be fed in instead of plasmid library.",
    )
    parser.add_argument(
        "--edit-start-pos",
        "-s",
        help="0-based start posiiton (inclusive) of edit relative to the start of guide spacer.",
    )
    parser.add_argument(
        "--edit-end-pos",
        "-e",
        help="0-based end position (exclusive) of edit relative to the start of guide spacer.",
    )
    parser.add_argument(
        "--jaccard-threshold",
        "-j",
        help="Jaccard Index threshold when the alleles are mapped to the most similar alleles. In each filtering step, allele counts of filtered out alleles will be mapped to the most similar allele only if they have Jaccard Index of shared edit higher than this threshold.",
        default=0.5,
    )
    parser.add_argument(
        "--filter_window",
        "-w",
        help="Only consider edit within window provided by (edit-start-pos, edit-end-pos). If this flag is not provided, `--edit-start-pos` and `--edit-end-pos` flags are ignored.",
        action="store_true",
    )
    parser.add_argument(
        "--filter_target_basechange",
        "-b",
        help="Only consider target edit (stored in bdata.uns['target_base_change'])",
        action="store_true",
    )
    parser.add_argument(
        "--translate", "-t", help="Translate alleles", action="store_true"
    )
    parser.add_argument(
        "--filter-allele-proportion",
        "-ap",
        type=float,
        default=None,
        help="If provided, alleles that exceed `filter_allele_proportion` in `filter-sample-proportion` will be retained.",
    )
    parser.add_argument(
        "--filter-sample-proportion",
        "-ap",
        type=float,
        default=0.2,
        help="If `filter_allele_proportion` is provided, alleles that exceed `filter_allele_proportion` in `filter-sample-proportion` will be retained.",
    )
    parser.add_argument(
        "--jaccard-threshold",
        "-j",
        help="Jaccard Index threshold when the alleles are mapped to the most similar alleles. In each filtering step, allele counts of filtered out alleles will be mapped to the most similar allele only if they have Jaccard Index of shared edit higher than this threshold.",
        default=0.5,
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    bdata = be.read_h5ad(args.h5ad_path)
    allele_df_keys = ["allele_counts"]
    if args.plasmid_path is not None:
        plasmid_adata = be.read_h5ad(args.plasmid_path)
        q_val_each, bdata.uns["sig_allele_counts"] = be.filter_alleles(
            bdata, plasmid_adata, filter_each_sample=True, run_parallel=True
        )
        allele_df_keys.append("sig_allele_counts")
    bdata.uns[f"{allele_df_keys[-1]}_spacer"] = bdata.filter_allele_counts_by_pos(
        rel_pos_start=0,
        rel_pos_end=20,
        rel_pos_is_reporter=False,
        map_to_filtered=True,
        allele_uns_key=allele_df_keys[-1],
        jaccard_threshold=args.jaccard_threshold,
    )
    allele_df_keys.append(f"{allele_df_keys[-1]}_spacer")
    if args.filter_window:
        filtered_key = f"{allele_df_keys[-1]}_{args.edit_start_pos}_{args.edit_end_pos}"
        bdata.uns[filtered_key] = bdata.filter_allele_counts_by_pos(
            rel_pos_start=args.edit_start_pos,
            rel_pos_end=args.edit_end_pos,
            rel_pos_is_reporter=False,
            map_to_filtered=True,
            allele_uns_key=allele_df_keys[-1],
            jaccard_threshold=args.jaccard_threshold,
        )

        allele_df_keys.append(filtered_key)
    if args.filter_target_basechange:
        bdata.uns[
            f"{allele_df_keys[-1]}_{bdata.uns['target_base_change']}"
        ] = bdata.filter_allele_counts_by_base(
            bdata.base_edited_from,
            bdata.base_edited_to,
            map_to_filtered=True,
            allele_uns_key=allele_df_keys[-1],
            jaccard_threshold=args.jaccard_threshold,
        )
        allele_df_keys.append(f"{allele_df_keys[-1]}_{bdata.uns['target_base_change']}")
    if args.translate:
        bdata.uns[f"{allele_df_keys[-1]}_translated"] = be.translate_allele_df(
            bdata.uns[f"{allele_df_keys[-1]}_translated"]
        ).rename(columns={"allele": "aa_allele"})
    if args.filter_allele_proportion is not None:
        filtered_key = f"{allele_df_keys[-1]}_prop{args.filter_allele_proportion}_{args.filter_sample_proportion}"
        bdata.uns[filtered_key] = be.an.filter_alleles.filter_allele_prop(
            bdata,
            allele_df_keys[-1],
            allele_prop_thres=args.filter_allele_proportion,
            sample_prop_thres=args.filter_sample_proportion,
            map_to_filtered=True,
            retain_max=True,
            allele_col=bdata.uns[allele_df_keys[-1]].columns[1],
            distribute=True,
            jaccard_threshold=args.jaccard_threshold,
        )
        allele_df_keys.append(filtered_key)

    if args.plot:
        fig, ax = plt.subplots(
            2, len(allele_df_keys), figsize=(6, 3 * len(allele_df_keys))
        )
        for i, key in enumerate(allele_df_keys):
            be.pl.allele_stats.plot_n_alleles_per_guide(
                bdata, key, bdata.uns[key].columns[1], ax[i, 0]
            )
            be.pl.allele_stats.plot_n_guides_per_edit(
                bdata, key, bdata.uns[key].columns[1], ax[i, 1]
            )
        plt.tight_layout()
        plt.savefig(
            f"{args.output_prefix}.filtered_allele_stats.pdf", bbox_inches="tight"
        )
    with open(f"{args.output_prefix}.filter_log.txt", "w") as out_log:
        for key in allele_df_keys:
            out_log.write(f"{key}\t{len(bdata.uns[key])}\n")
