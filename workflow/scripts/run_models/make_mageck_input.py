import os
import argparse
import patsy
import pandas as pd
import bean as be
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("bdata_path", type=str, help="Path of an ReporterScreen object")
parser.add_argument("--prefix", "-p", default="", help="prefix to save the result")
parser.add_argument(
    "--dmatrix_topbot", "-t", default=True, help="Write top/bot design matrix"
)
parser.add_argument(
    "--dmatrix_sort", "-s", default=True, help="Write sort design matrix"
)
parser.add_argument(
    "--sample_mask_column",
    "-m",
    help="column in ReporterScreen.samples specifying which samples will be used. Sample is used if 1 and not used if 0.",
)
parser.add_argument(
    "--reps",
    type=str,
    nargs="+",
    default=[],
    help="Subset these replicates for the run.",
)
parser.add_argument(
    "--use_bcmatch", action="store_true", help="Use barcode matched reads to train"
)
parser.add_argument(
    "--target_col",
    "-c",
    type=str,
    default="target",
    help="Column name of bdata.guides that shows the target element of each guide",
)
args = parser.parse_args()

matrix_prefix = (
    args.prefix + "/" + os.path.basename(args.bdata_path).rsplit(".h5ad", 1)[0]
)
if args.reps:
    args.prefix = (
        args.prefix
        + ".rep"
        + "_".join([rep.replace("_", "") for rep in args.reps])
    )

bdata = be.read_h5ad(args.bdata_path)
if args.reps:
    print(f"Subsetting the specified replicates {args.reps}")
    bdata = bdata[:, bdata.samples.rep.isin(args.reps)]
if args.use_bcmatch:
    Xdf = bdata.to_mageck_input(target_column=args.target_col)
    bcmatch_df = bdata.to_mageck_input(
        count_layer="X_bcmatch", target_column=args.target_col, sample_prefix="bcmatch_"
    )
    Xdf.merge(bcmatch_df, on=["sgRNA", "gene"]).to_csv(
        f"{matrix_prefix}.bcmatch.mageck_input.txt", sep="\t", index=False
    )
else:
    bdata.to_mageck_input(
        out_path="{}.mageck_input.txt".format(matrix_prefix),
        target_column=args.target_col,
    )
edit_rates = bdata.guides["edit_rate"].fillna(0.5)
edit_rates.to_csv("{}.sgrna_eff.txt".format(matrix_prefix), header=False, sep="\t")
(edit_rates * 1.25 - 1).to_csv(
    "{}.mageck_sgrna_eff.txt".format(matrix_prefix), header=False, sep="\t"
)
if not args.sample_mask_column is None:
    bdata_use = bdata[:, bdata.samples[args.sample_mask_column] == 1]
else:
    bdata_use = bdata

if args.dmatrix_topbot:
    if not "sort" in bdata_use.samples.columns:
        bdata_use.samples["sort"] = bdata_use.samples.index.map(
            lambda s: s.rsplit("_")[-1]
        )
    bdata_topbot = bdata_use[:, bdata_use.samples.sort.isin(["top", "bot"])]
    dm = patsy.dmatrix("sort", bdata_topbot.samples, return_type="dataframe")
    dm = dm.set_index(bdata_topbot.samples.index).reset_index()
    dm = dm.rename(
        columns={
            "index": "Samples",
            "Intercept": "baseline",
            "sort[T.top]": "top_vs_bot",
        }
    )
    with warnings.catch_warnings():
        dm.iloc[:, 1:] = dm.iloc[:, 1:].astype(int)
    if args.use_bcmatch:
        dm_bcmatch = dm.set_index("Samples").rename("bcmatch_{}".format).reset_index()
        dm = pd.concat((dm, dm_bcmatch), ignore_index=True)
    dm.to_csv("{}.mageck_dm_topbot.txt".format(matrix_prefix), sep="\t", index=False)

if args.dmatrix_sort:
    if not "sort" in bdata_use.samples.columns:
        bdata_use.samples["sort"] = bdata_use.samples.index.map(
            lambda s: s.rsplit("_")[-1]
        )
    bdata_sort = bdata_use[:, bdata_use.samples.sort != "bulk"]
    bdata_sort.samples["sort_num"] = 0
    bdata_sort.samples.loc[bdata_sort.samples.sort == "low", "sort_num"] = 1
    bdata_sort.samples.loc[bdata_sort.samples.sort == "high", "sort_num"] = 3
    bdata_sort.samples.loc[bdata_sort.samples.sort == "top", "sort_num"] = 4
    dm = patsy.dmatrix("sort_num", bdata_sort.samples, return_type="dataframe")
    dm = dm.set_index(bdata_sort.samples.index).reset_index()
    dm = dm.rename(
        columns={
            "index": "Samples",
            "Intercept": "baseline",
            "sort[T.top]": "top_vs_bot",
        }
    )
    dm.iloc[:, 1:] = dm.iloc[:, 1:].astype(int)
    if args.use_bcmatch:
        dm_bcmatch = dm.set_index("Samples").rename("bcmatch_{}".format).reset_index()
        dm = pd.concat((dm, dm_bcmatch), ignore_index=True)
    dm.to_csv("{}.mageck_dm_sort.txt".format(matrix_prefix), sep="\t", index=False)
