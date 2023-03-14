import os
from copy import deepcopy
import numpy as np
import torch
import pyro
import pandas as pd
import pickle as pkl
import bean as be
from utils import (
    _get_guide_target_info,
    run_inference,
    run_inference_custom_guide,
    run_MCMC,
    write_result_table,
    get_parser_args,
)

os.environ["OPENBLAS_NUM_THREADS"] = "1"
pyro.set_rng_seed(101)


def main(args):
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

    bdata = be.read_h5ad(args.bdata_path)

    if not args.accessibility_bw is None and (not "chr" in bdata.guides.columns):
        updated_guide_info = pd.read_csv(
            "/data/pinello/PROJECTS/2021_08_ANBE/data/gRNA_info/LDLvar_gRNA_beret_accessibility.csv",
            index_col=0,
        )
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
    guide_index = bdata.guides.index
    if args.reps:
        print(f"Subsetting the specified replicates {args.reps}")
        bdata = bdata[:, bdata.samples.rep.isin(args.reps)]

    bdata_fullsort = bdata[:, bdata.samples.sorting_scheme == "sort"]
    bdata_topbot = bdata[:, bdata.samples.sorting_scheme == "topbot"]
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

    if bdata_topbot.shape[1] != 0:
        if not args.repguide_mask is None:
            assert args.repguide_mask in bdata_topbot.uns.keys()
        ndata_topbot = load_data_normal.NData(
            bdata_topbot,
            5,
            "topbot",
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
        ndata_fullsort = load_data_normal.NData(
            bdata_fullsort,
            5,
            "bins",
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
        negctrl_idx = np.where(bdata.guides.Group.str.lower() == "negctrl")[0]

    if args.accessibility_bw or args.accessibility_col:
        bdata.guides["accessibility"] = ndatas[0].guide_accessibility
        target_info = _get_guide_target_info(bdata, guide_acc_col="accessibility")
    else:
        target_info = _get_guide_target_info(bdata)
    if ndatas[0].guide_accessibility is None:
        guide_acc = None
    else:
        guide_acc = ndatas[0].guide_accessibility.detach().numpy()
    outfile_path = "{}.model{}.result.pkl".format(args.prefix, args.model_id)
    if not os.path.exists(outfile_path):
        del bdata
        del bdata_fullsort
        del bdata_topbot
        use_bcmatch = args.use_bcmatch
        model_dict = {
            "A": lambda data: m.model_A_both(data, use_bcmatch=use_bcmatch),
            "B0": lambda data: m.model_B0_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B_sp": lambda data: m.model_B_samplepi_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B_spm": lambda data: m.model_B_samplepi_multinom_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B": lambda data: m.model_B_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B1": lambda data: m.model_B1_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B2": lambda data: m.model_B2_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B3": lambda data: m.model_B3_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B4": lambda data: m.model_B4_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            # "A_sl": lambda data: msl.model_A_both(data, use_bcmatch = use_bcmatch),
            # "B0_sl": lambda data: msl.model_B0_both(data, alpha_prior=10, use_bcmatch = use_bcmatch),
            # "B_sl": lambda data: msl.model_B_both(data, alpha_prior=10, use_bcmatch = use_bcmatch),
            # "B2_sl": lambda data: msl.model_B2_both(data, alpha_prior=10, use_bcmatch = use_bcmatch),
            "ANB": lambda data: mnb.model_A_NB_both(data, use_bcmatch=use_bcmatch),
            "B0NB": lambda data: mnb.model_B0_NB_both(data, use_bcmatch=use_bcmatch),
            "BNB": lambda data: mnb.model_B_NB_both(data, use_bcmatch=use_bcmatch),
            # "BNB_sp": lambda data: mnb.model_B_NB_sample_pis(data, use_bcmatch = use_bcmatch),
            "B2NB": lambda data: mnb.model_B2_NB_both(data, use_bcmatch=use_bcmatch),
            "B5": lambda data: m.model_B5_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B10": lambda data: m.model_B10_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B14": lambda data: m.model_B14_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B14_0": lambda data: m.model_B14_both_0(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B14_2": lambda data: m.model_B14_both_2(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B15": lambda data: m.model_B15_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B16": lambda data: m.model_B16_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B16_2": lambda data: m.model_B16_both_2(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B16_0": lambda data: m.model_B16_both_0(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B16_1": lambda data: m.model_B16_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B17": lambda data: m.model_B17_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B18": lambda data: m.model_B18_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B19": lambda data: m.model_B19_both(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            "B18_2": lambda data: m.model_B18_both_2(
                data, alpha_prior=10, use_bcmatch=use_bcmatch
            ),
            # "B2NB_sp": lambda data: mnb.model_B2_NB_sample_pis(data, use_bcmatch = use_bcmatch),
        }

        model = model_dict[args.model_id]
        guide_dict = {
            "B14": m.guide_B14_both,
            "B14_2": m.guide_B14_both_2,
            "B14_0": m.guide_B14_both,
            "B16": m.guide_B14_both,
            "B16_1": m.guide_B14_both,
            "B16_2": m.guide_B14_both_2,
            "B16_0": m.guide_B14_both,
            "B15": m.guide_B_both,
            "B18": m.guide_B14_both,
            "B18_2": m.guide_B14_both_2,
            "B19": m.guide_B14_both,
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

        print("Done running inference. Writing result at {}...".format(outfile_path))
        with open(
            "{}.model{}.result.pkl".format(args.prefix, args.model_id), "wb"
        ) as handle:
            pkl.dump(param_history_dict, handle)

        write_result_table(
            target_info,
            param_history_dict,
            prefix=f"{args.prefix}.model{args.model_id}.",
            write_fitted_eff=(args.model_id != "A" and not args.model_id.endswith("0")),
            guide_index=guide_index,
            guide_acc=guide_acc,
        )

    else:
        print(f"Reading saved result at {outfile_path}:")
        with open(outfile_path, "rb") as f:
            param_history_dict = pkl.load(f)
        write_result_table(
            target_info,
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
