import os
from copy import deepcopy
import numpy as np
import torch
import pyro
import pandas as pd
import pickle as pkl
import bean as be
from utils import (
    run_inference,
    run_inference_custom_guide,
    run_MCMC,
    write_result_table,
    get_parser_args,
    annotate_edit,
)

os.environ["OPENBLAS_NUM_THREADS"] = "1"
pyro.set_rng_seed(101)


def main(args):
    bdata = be.read_h5ad(args.bdata_path)
    guide_index = bdata.guides.index

    if args.cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)
    else:
        torch.set_default_tensor_type(torch.FloatTensor)
    if args.prefix == "":
        args.prefix = "."
    os.makedirs(args.prefix, exist_ok=True)
    args.prefix = (
        args.prefix + "/" + os.path.basename(args.bdata_path).rsplit(".h5ad", 1)[0]
    )
    if args.reps:
        args.prefix = (
            args.prefix
            + f".rep"
            + "_".join([rep.replace("_", "") for rep in args.reps])
        )
    outfile_path = "{}.model{}.result.pkl".format(args.prefix, args.model_id)
    if not args.accessibility_bw is None and (
        not args.accessibility_col in bdata.guides.columns
    ):
        updated_guide_info = pd.read_csv(
            "/data/pinello/PROJECTS/2021_08_ANBE/data/gRNA_info/LDLRCDS_gRNA_beret_accessibility.csv",
            index_col=0,
        )  # TODO: make this file
        bdata.guides = bdata.guides.join(
            updated_guide_info[["genomic_pos", "chr", "parsed_snp"]]
        )
        print("Writing bdata with updated target info")
        bdata.write(args.bdata_path)

    if not "edit_rate" in bdata.guides.columns:
        bdata.guides["edit_rate"] = bdata.get_guide_edit_rate(return_result=True)
        bdata.write(args.bdata_path)
    print("Done loading data. Preprocessing...")
    if args.model_id in ["B16_0", "B0"]:
        # For models that doesn't use observed editing rates to fit pi_a0
        bdata.layers["edits"] = (bdata.layers["X_bcmatch"] / 2.0).astype(int)

    if args.reps:
        print(f"Subsetting the specified replicates {args.reps}")
        bdata = bdata[:, bdata.samples.rep.isin(args.reps)]
    bdata_fullsort = bdata[:, bdata.samples.sorting_scheme == "sort"]
    bdata_topbot = bdata[:, bdata.samples.sorting_scheme == "topbot"]
    bdata_bulk = bdata[:, bdata.samples.sort == "bulk"]
    bdata_bulk.remove_zero_allele_counts()

    alleles_to_select = (
        bdata_bulk.uns[args.allele_df_key]
        .set_index(bdata_bulk.uns[args.allele_df_key].columns[:2].tolist())
        .index
    )
    ndatas = []
    bdata.samples["rep"] = bdata.samples["rep"].astype("category")
    total_reps = bdata.samples["rep"].unique()
    device = None

    if args.mcmc:
        device = "cuda:0"
    dispersion_path = None
    if "NB" in args.model_id and args.fit_dispersion:
        tmpfile_prefix = args.bdata_path.rsplit(".", 1)[0]
        bdata_X_path = tmpfile_prefix + ".X.txt"
        bdata_condit_path = tmpfile_prefix + ".samples.txt"
        dispersion_path = args.bdata_path.rsplit(".", 1)[0] + ".DESeq2_dispersion.csv"
        if not os.path.exists(dispersion_path):
            np.savetxt(bdata_X_path, bdata.X, delimiter=",")
            bdata.samples.to_csv(bdata_condit_path)
            os.system(
                "Rscript ../get_DESeq_dispersion.R {} {} {}".format(
                    bdata_X_path, bdata_condit_path, dispersion_path
                )
            )
            print("Fitted dispersion saved at {}.".format(dispersion_path))
        else:
            print("Using already fitted dispersion at {}.".format(dispersion_path))

    edit_index = load_data_normal.get_edit_to_index_dict(
        bdata_bulk.uns[args.allele_df_key].aa_allele
    )

    edit_index_df = pd.DataFrame(pd.Series(edit_index))
    edit_index_df.to_csv("{}.model{}.edit_index.csv".format(args.prefix, args.model_id))
    if bdata_topbot.shape[1] != 0:
        if not args.repguide_mask is None:
            assert args.repguide_mask in bdata_topbot.uns.keys()
        ndata_topbot = load_data_normal.MData(
            bdata_topbot,
            5,
            "topbot",
            allele_df_key=args.allele_df_key,
            edit_index=edit_index,
            dispersion_path=dispersion_path,
            device=device,
            repguide_mask=args.repguide_mask,
            sample_mask_column=args.sample_mask_column,
            fit_a0=~args.raw_alpha,
            fit_a0_lower_quantile=args.a0_lower_quantile,
            accessibility_bw_path=args.accessibility_bw,
            accessibility_col=args.accessibility_col,
            shrink_alpha=args.shrink_alpha,
            impute_pi_popt=False,
            use_const_pi=args.model_id in ["B1", "B16_1"],
            alleles_to_select=alleles_to_select,
        )
        if args.mcmc:
            print(ndata_topbot.X.get_device())
        ndata_topbot.rep_index = np.where(
            total_reps.isin(bdata_topbot.samples["rep"].unique())
        )[0]
        ndatas.append(ndata_topbot)

    if bdata_fullsort.shape[1] != 0:
        if not args.repguide_mask is None:
            assert args.repguide_mask in bdata_fullsort.uns.keys()
        ndata_fullsort = load_data_normal.MData(
            bdata_fullsort,
            5,
            "bins",
            allele_df_key=args.allele_df_key,
            edit_index=edit_index,
            dispersion_path=dispersion_path,
            device=device,
            repguide_mask=args.repguide_mask,
            sample_mask_column=args.sample_mask_column,
            fit_a0=~args.raw_alpha,
            fit_a0_lower_quantile=args.a0_lower_quantile,
            accessibility_bw_path=args.accessibility_bw,
            accessibility_col=args.accessibility_col,
            shrink_alpha=args.shrink_alpha,
            impute_pi_popt=False,
            use_const_pi=args.model_id in ["B1", "B16_1"],
            alleles_to_select=alleles_to_select,
        )
        if args.mcmc:
            print(ndata_fullsort.X.get_device())  # would throw error if on cpu
        ndata_fullsort.rep_index = np.where(
            total_reps.isin(bdata_fullsort.samples["rep"].unique())
        )[0]
        ndatas.append(ndata_fullsort)

    for ndata in ndatas:
        ndata.n_total_reps = len(total_reps)
    if args.fit_negctrl:
        negctrl_idx = np.where(bdata.guides.index.map(lambda s: "CONTROL" in s))[0]
    del bdata
    del bdata_fullsort
    del bdata_topbot
    del bdata_bulk
    use_bcmatch = args.use_bcmatch
    model_dict = {
        # "A": lambda data: m.model_A_both(data, use_bcmatch=use_bcmatch),
        "C14": lambda data: m.model_C14_both(
            data, alpha_prior=10, use_bcmatch=use_bcmatch
        ),
        # "C14_0": lambda data: m.model_C14_both_0(
        #     data, alpha_prior=10, use_bcmatch=use_bcmatch
        # ),
        # "C14_2": lambda data: m.model_C14_both_2(
        #     data, alpha_prior=10, use_bcmatch=use_bcmatch
        # ),
        "C16": lambda data: m.model_C16_both(
            data, alpha_prior=10, use_bcmatch=use_bcmatch
        ),
        # "C16_2": lambda data: m.model_C16_both_2(
        #     data, alpha_prior=10, use_bcmatch=use_bcmatch
        # ),
        # "C16_0": lambda data: m.model_C16_both_0(
        #     data, alpha_prior=10, use_bcmatch=use_bcmatch
        # ),
        # "C16_1": lambda data: m.model_C16_both(
        #     data, alpha_prior=10, use_bcmatch=use_bcmatch
        # ),
        "C17": lambda data: m.model_C17_both(
            data, alpha_prior=10, use_bcmatch=use_bcmatch
        ),
        "C18": lambda data: m.model_C18_both(
            data, alpha_prior=10, use_bcmatch=use_bcmatch
        ),
        "C19": lambda data: m.model_C19_both(
            data, alpha_prior=10, use_bcmatch=use_bcmatch
        ),
        # "C18_2": lambda data: m.model_C18_both_2(
        #     data, alpha_prior=10, use_bcmatch=use_bcmatch
        # ),
    }
    if ndatas[0].guide_accessibility is None:
        guide_acc = None
    else:
        guide_acc = ndatas[0].guide_accessibility.detach().numpy()
    if not os.path.exists(outfile_path):
        print(f"Running model as {outfile_path} not found:")
        model = model_dict[args.model_id]
        guide_dict = {
            "C14": m.guide_C14_both,
            # "C14_2": m.guide_C14_both_2,
            # "C14_0": m.guide_C14_both,
            "C16": m.guide_C14_both,
            # "C16_1": m.guide_C14_both,
            # "C16_2": m.guide_C14_both_2,
            # "C16_0": m.guide_C14_both,
            "C18": m.guide_C14_both,
            # "C18_2": m.guide_C14_both_2,
            "C19": m.guide_C14_both,
        }
        print("Running inference for model {}...".format(args.model_id))
        if args.mcmc:
            param_history_dict = run_MCMC(model, ndatas)
        else:
            if args.fit_negctrl:
                negctrl_model = m.model_A0_both
                negctrl_guide = m.guide_A0_both
                ndatas_negctrl = [ndata[negctrl_idx] for ndata in ndatas]
                param_dict_negctrl = deepcopy(
                    run_inference_custom_guide(
                        negctrl_model, negctrl_guide, ndatas_negctrl
                    )
                )
            if args.model_id in guide_dict.keys():
                param_history_dict = run_inference_custom_guide(
                    model, guide_dict[args.model_id], ndatas
                )
            else:
                print(
                    f"No guide specified for model {args.model_id}. Using AutoNormal guide"
                )
                param_history_dict = run_inference(
                    model,
                    ndatas,
                    run_id=f"{os.path.basename(args.bdata_path)}_{args.model_id}",
                )
            if args.fit_negctrl:
                param_history_dict["negctrl"] = param_dict_negctrl
        param_history_dict["edit_index"] = edit_index_df.reset_index()

        print("Done running inference. Writing result at {}...".format(outfile_path))
        with open(
            "{}.model{}.result.pkl".format(args.prefix, args.model_id), "wb"
        ) as handle:
            pkl.dump(param_history_dict, handle)

    else:
        print(f"Reading saved result at {outfile_path}:")
        with open(outfile_path, "rb") as f:
            param_history_dict = pkl.load(f)
        try:
            edit_index_df = (
                param_history_dict["edit_index"]
                .reset_index()
                .rename(columns={"index": "edit"})
            )
            edit_index_df = annotate_edit(edit_index_df)
        except Exception as e:
            print(e)
            edit_index_df = param_history_dict["edit_index"].reset_index()
    write_result_table(
        edit_index_df,  # TODO: don't need reset index later on
        param_history_dict,
        prefix=f"{args.prefix}.model{args.model_id}.",
        write_fitted_eff=(args.model_id != "A" and not args.model_id.endswith("0")),
        guide_index=guide_index,
        guide_acc=guide_acc,
    )
    print("Done!")


if __name__ == "__main__":
    args = get_parser_args()
    main(args)
