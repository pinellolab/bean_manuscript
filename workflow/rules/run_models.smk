
rule run_bean_var_normal:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated{suffix}.h5ad"
    output:
        output="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated{suffix}/bean_element_result.Normal.csv"
    run:
        shell("bean-run variant {input.input_h5ad} --uniform-edit -o results/model_runs/bean/ --ignore-bcmatch")

rule run_bean_var_mixnormal:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated{suffix}.h5ad"
    output:
        output="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated{suffix}/bean_element_result.MixtureNormal.csv"
    run:
        shell("bean-run variant {input.input_h5ad} -o results/model_runs/bean/ --ignore-bcmatch")

rule run_bean_var_mixnormal_acc:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated{suffix}.h5ad"
    output:
        output="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated{suffix}/bean_element_result.MixtureNormal+Acc.csv"
    run:
        shell("bean-run variant {input.input_h5ad} --scale-by-acc --acc-bw-path resources/accessibility/ENCFF262URW.hg19.bw -o results/model_runs/bean/ --ignore-bcmatch")

rule run_bean_var_normal_negctrl:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated{suffix}.h5ad"
    output:
        output="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLvar_annotated{suffix,.*}/bean_element_result.Normal.csv"
    run:
        shell("bean-run variant {input.input_h5ad} --uniform-edit -o results/model_runs/bean_negctrl/ --fit-negctrl  ")

rule run_bean_var_mixnormal_negctrl:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated{suffix}.h5ad"
    output:
        output="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLvar_annotated{suffix,.*}/bean_element_result.MixtureNormal.csv"
    run:
        shell("bean-run variant {input.input_h5ad} -o results/model_runs/bean_negctrl/ --fit-negctrl  ")

rule run_bean_var_mixnormal_acc_negctrl:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated{suffix}.h5ad"
    output:
        output="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLvar_annotated{suffix,.*}/bean_element_result.MixtureNormal+Acc.csv"
    run:
        shell("bean-run variant {input.input_h5ad} --scale-by-acc --acc-bw-path resources/accessibility/ENCFF262URW.hg19.bw -o results/model_runs/bean_negctrl/ --fit-negctrl  ")

rule run_mageck_var:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated_complete.h5ad"
    output:
        mle_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/sort.gene_summary.txt",
        mle_var_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/sort_var.gene_summary.txt",
        mle_em_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/EM/sort.gene_summary.txt",
        mle_var_em_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/EM/sort_var.gene_summary.txt",
        mle_emf_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/EMf/sort.gene_summary.txt",
        mle_var_emf_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/EMf/sort_var.gene_summary.txt",
        rra_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/rra_top.gene_summary.txt",
        rra_bot_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/rra_bot.gene_summary.txt",
    run:
        shell("sh scripts/run_models/run_mageck_var.sh {input.input_h5ad} results/model_runs/mageck")

rule run_mageck_var_negctrl:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated_complete.h5ad",
        negctrl_guides="resources/gRNA_info/LDLvar_negctrl_gRNA.txt"
    output:
        mle_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/sort.gene_summary.txt",
        mle_var_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/sort_var.gene_summary.txt",
        mle_em_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/EM/sort.gene_summary.txt",
        mle_var_em_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/EM/sort_var.gene_summary.txt",
        mle_emf_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/EMf/sort.gene_summary.txt",
        mle_var_emf_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/EMf/sort_var.gene_summary.txt",
        rra_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/rra_top.gene_summary.txt",
        rra_bot_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/rra_bot.gene_summary.txt",
    run:
        shell("sh scripts/run_models/run_mageck_var.sh {input.input_h5ad} results/model_runs/mageck_negctrl --control-sgrna {input.negctrl_guides}")

rule run_cb2:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated_complete.h5ad"
    output:
        cb2_out="results/model_runs/CB2/CB2_run_result.bean_count_LDLvar_annotated_complete/CB2_gene.csv",
    run:
        shell("conda run -n anbe_benchmark sh scripts/run_models/run_CB2_var.sh {input.input_h5ad}")
        

rule evaluate_varscreen:
    input:
        mle_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/sort.gene_summary.txt",
        mle_var_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/sort_var.gene_summary.txt",
        mle_em_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/EM/sort.gene_summary.txt",
        mle_var_em_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/EM/sort_var.gene_summary.txt",
        mle_emf_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/EMf/sort.gene_summary.txt",
        mle_var_emf_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/EMf/sort_var.gene_summary.txt",
        rra_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/rra_top.gene_summary.txt",
        rra_bot_out="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/rra_bot.gene_summary.txt",
        normal_res="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated_complete/bean_element_result.Normal.csv",
        mixnormal_res="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated_complete/bean_element_result.MixtureNormal.csv",
        mixnormal_acc_res="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated_complete/bean_element_result.MixtureNormal+Acc.csv",
        cb2_res="results/model_runs/CB2/CB2_run_result.bean_count_LDLvar_annotated_complete/CB2_gene.csv",
    params:
        bean_prefix="results/model_runs/bean/bean_count_LDLvar_annotated_complete",
        mageck_prefix="results/model_runs/mageck/bean_count_LDLvar_annotated_complete/",
    output:
        "results/model_runs/bean_count_LDLvar_annotated_complete/all_scores.csv",
    run:
        shell("python scripts/evaluate_model_runs/get_performance_varscreen.py bean_count_LDLvar_annotated_complete")
        shell("python scripts/evaluate_model_runs/get_performance_varscreen.py bean_count_LDLvar_annotated_complete --control strongest_splicing")

rule evaluate_varscreen_negctrl:
    input:
        mle_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/sort.gene_summary.txt",
        mle_var_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/sort_var.gene_summary.txt",
        mle_em_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/EM/sort.gene_summary.txt",
        mle_var_em_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/EM/sort_var.gene_summary.txt",
        mle_emf_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/EMf/sort.gene_summary.txt",
        mle_var_emf_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/EMf/sort_var.gene_summary.txt",
        rra_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/rra_top.gene_summary.txt",
        rra_bot_out="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/rra_bot.gene_summary.txt",
        normal_res="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLvar_annotated_complete/bean_element_result.Normal.csv",
        mixnormal_res="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLvar_annotated_complete/bean_element_result.MixtureNormal.csv",
        mixnormal_acc_res="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLvar_annotated_complete/bean_element_result.MixtureNormal+Acc.csv",
        
    params:
        bean_prefix="results/model_runs/bean_negctrl/bean_count_LDLvar_annotated_complete",
        mageck_prefix="results/model_runs/mageck_negctrl/bean_count_LDLvar_annotated_complete/",
    output:
        "results/model_runs/bean_count_LDLvar_annotated_complete/all_scores_negctrl.csv",
    run:
        shell("python scripts/evaluate_model_runs/get_performance_varscreen.py bean_count_LDLvar_annotated_complete --result-suffix _negctrl")
        shell("python scripts/evaluate_model_runs/get_performance_varscreen.py bean_count_LDLvar_annotated_complete --control strongest_splicing --result-suffix _negctrl")


rule run_2reps_varscreen:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated_complete.h5ad"
    output:
        "results/model_runs/bean_count_LDLvar_annotated_complete_rep14_rep15/all_scores.csv"
    run:
        shell("python scripts/run_models/run_models_on_2_reps.py {input.input_h5ad}")


rule run_bean_tiling_negctrl:
    input:
        input_h5ad="results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_annotated{cutoff_suffix}.h5ad",
        splice_site="resources/LDLR/LDLR_ABE_splice_targets.csv"
    params:
        trailing_args= lambda wildcards: "--ignore-bcmatch" if "complete" in wildcards.cutoff_suffix else ""
    output:
        res="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLRCDS_annotated{cutoff_suffix}/bean_element_result.MultiMixtureNormal.csv",
    run:
        shell("bean-run tiling {input.input_h5ad} -o results/model_runs/bean_negctrl --allele-df-key sig_allele_counts_spacer_0_19_A.G_translated_prop0.1_0.3 --splice-site-path {input.splice_site} --control-guide-tag ABE_CONTROL --fit-negctrl --negctrl-col Region --negctrl-col-value 'ABE control' --cuda {params.trailing_args}")

rule run_bean_tiling_negctrl_norm:
    input:
        input_h5ad="results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_annotated{cutoff_suffix}.h5ad",
    params:
        trailing_args= lambda wildcards: "--ignore-bcmatch" if "complete" in wildcards.cutoff_suffix else ""
    output:
        normal_res_all="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLRCDS_annotated{cutoff_suffix}/bean_element_result.Normal_allEdited.csv",
        normal_res_pred="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLRCDS_annotated{cutoff_suffix}/bean_element_result.Normal_behive.csv",
    run:
        shell("bean-run variant {input.input_h5ad} --uniform-edit -o results/model_runs/bean_negctrl --control-guide-tag ABE_CONTROL --fit-negctrl --negctrl-col Region --negctrl-col-value 'ABE control' --target-col target_allEdited --result-suffix _allEdited --cuda {params.trailing_args}")
        shell("bean-run variant {input.input_h5ad} --uniform-edit -o results/model_runs/bean_negctrl --control-guide-tag ABE_CONTROL --fit-negctrl --negctrl-col Region --negctrl-col-value 'ABE control' --target-col target_behive --result-suffix _behive --cuda {params.trailing_args}")

rule run_bean_tiling_acc_negctrl:
    input:
        input_h5ad="results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_annotated{cutoff_suffix}.h5ad",
        splice_site="resources/LDLR/LDLR_ABE_splice_targets.csv"
    params:
        trailing_args= lambda wildcards: "--ignore-bcmatch" if "complete" in wildcards.cutoff_suffix else ""
    output:
        res="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLRCDS_annotated{cutoff_suffix}/bean_element_result.MultiMixtureNormal+Acc.csv",
    run:
        shell("bean-run tiling {input.input_h5ad} --scale-by-acc --acc-bw-path resources/accessibility/ENCFF262URW.hg19.bw -o results/model_runs/bean_negctrl --allele-df-key sig_allele_counts_spacer_0_19_A.G_translated_prop0.1_0.3 --splice-site-path {input.splice_site} --control-guide-tag ABE_CONTROL --fit-negctrl --negctrl-col Region --negctrl-col-value 'ABE control' --cuda {params.trailing_args}")

rule run_bean_tiling_cbe_negctrl:
    input:
        input_h5ad="results/filtered_annotated/LDLRCDS_CBE_{cas_enyme}/bean_count_LDLRCDS_CBE_{cas_enyme}_annotated{cutoff_suffix}.h5ad",
        splice_site="resources/LDLR/LDLR_CBE_splice_targets.csv"
    output:
        res="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLRCDS_CBE_{cas_enyme}_annotated/bean_element_result.MultiMixtureNormal.csv",
    run:
        shell("bean-run tiling {input.input_h5ad} -o results/model_runs/bean_negctrl --allele-df-key sig_allele_counts_spacer_3_8_C.T_translated_prop0.1_0.3 --splice-site-path {input.splice_site} --control-guide-tag CBE_CONTROL --cuda --fit-negctrl --negctrl-col Region --negctrl-col-value 'CBE control'")

rule run_bean_tiling_acc_cbe_negctrl:
    input:
        input_h5ad="results/filtered_annotated/LDLRCDS_CBE_{cas_enyme}/bean_count_LDLRCDS_CBE_{cas_enyme}_annotated{cutoff_suffix}.h5ad",
        splice_site="resources/LDLR/LDLR_CBE_splice_targets.csv"
    output:
        res="results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLRCDS_CBE_{cas_enyme}_annotated{cutoff_suffix}/bean_element_result.MultiMixtureNormal+Acc.csv",
    run:
        shell("bean-run tiling {input.input_h5ad} --scale-by-acc --acc-bw-path resources/accessibility/ENCFF262URW.hg19.bw -o results/model_runs/bean_negctrl --allele-df-key sig_allele_counts_spacer_3_8_C.T_translated_prop0.1_0.3 --splice-site-path {input.splice_site} --control-guide-tag CBE_CONTROL --cuda --fit-negctrl --negctrl-col Region --negctrl-col-value 'CBE control'")


rule run_mageck_tiling_negctrl:
    input:
        input_h5ad="results/filtered_annotated/{tiling_lib}/bean_count_{tiling_lib}_annotated_complete.h5ad"
    output:
        mle_out_all="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_allEdited/sort.gene_summary.txt",
        mle_out_pred="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_behive/sort.gene_summary.txt",
        mle_var_out_all="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_allEdited/sort_var.gene_summary.txt",
        mle_var_out_pred="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_behive/sort_var.gene_summary.txt",
        rra_out_all="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_allEdited/rra_top.gene_summary.txt",
        rra_out_pred="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_behive/rra_top.gene_summary.txt",
        rra_bot_out_all="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_allEdited/rra_bot.gene_summary.txt",
        rra_bot_out_pred="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_behive/rra_bot.gene_summary.txt",
    run:
        shell("sh scripts/run_models/run_mageck_tiling.sh {input.input_h5ad} results/model_runs/mageck_negctrl/ --negctrl")

rule run_cb2_tiling:
    input:
        input_h5ad="results/filtered_annotated/{tiling_lib}/bean_count_{tiling_lib}_annotated_complete.h5ad"
    output:
        cb2_out_all="results/model_runs/CB2/CB2_run_result.bean_count_{tiling_lib}_annotated_complete.target_allEdited/CB2_gene.csv",
        cb2_out_pred="results/model_runs/CB2/CB2_run_result.bean_count_{tiling_lib}_annotated_complete.target_behive/CB2_gene.csv",
    run:
        shell("conda run -n anbe_benchmark sh scripts/run_models/run_CB2_tiling.sh {input.input_h5ad}")

rule run_crisphieRmix_tiling_negctrl:
    input:
        input_h5ad="results/filtered_annotated/{tiling_lib}/bean_count_{tiling_lib}_annotated_complete.h5ad"
    params:
        control_label=lambda wildcards: "CBE control" if "CBE" in wildcards.tiling_lib else "ABE control"
    output:
        crisphiermix_out="results/model_runs/CRISPhieRmix_negctrl/CRISPhieRmix_run_result.bean_count_{tiling_lib}_annotated_complete.target_behive/CRISPhieRmix.csv",
        crisphiermix_out_all="results/model_runs/CRISPhieRmix_negctrl/CRISPhieRmix_run_result.bean_count_{tiling_lib}_annotated_complete.target_allEdited/CRISPhieRmix.csv",
    run:
        shell('mkdir -p results/model_runs/CRISPhieRmix_negctrl/CRISPhieRmix_run_result.bean_count_{wildcards.tiling_lib}_annotated_complete.target_behive/; mkdir -p results/model_runs/CRISPhieRmix_negctrl/CRISPhieRmix_run_result.bean_count_{wildcards.tiling_lib}_annotated_complete.target_allEdited/;  conda run -n anbe_benchmark Rscript scripts/run_models/run_CRISPhieRmix_LDLRCDS.R -i {input.input_h5ad} -o results/model_runs/CRISPhieRmix_negctrl/CRISPhieRmix_run_result.bean_count_{wildcards.tiling_lib}_annotated_complete')
        #shell("conda run -n anbe_benchmark sh scripts/run_models/run_CRISPhieRmix_tiling_negctrl.sh {input.input_h5ad} {params.control_label}")


rule collect_scores_tiling_negctrl:
    input:
        mle_out_all="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_allEdited/sort.gene_summary.txt",
        mle_out_pred="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_behive/sort.gene_summary.txt",
        mle_var_out_all="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_allEdited/sort_var.gene_summary.txt",
        mle_var_out_pred="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_behive/sort_var.gene_summary.txt",
        rra_out_all="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_allEdited/rra_top.gene_summary.txt",
        rra_out_pred="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_behive/rra_top.gene_summary.txt",
        rra_bot_out_all="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_allEdited/rra_bot.gene_summary.txt",
        rra_bot_out_pred="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete.no_bcmatch.target_behive/rra_bot.gene_summary.txt",
        normal_res_all="results/model_runs/bean_negctrl/bean_run_result.bean_count_{tiling_lib}_annotated_complete/bean_element_result.Normal_allEdited.csv",
        normal_res_pred="results/model_runs/bean_negctrl/bean_run_result.bean_count_{tiling_lib}_annotated_complete/bean_element_result.Normal_behive.csv",
        mixnormal_res="results/model_runs/bean_negctrl/bean_run_result.bean_count_{tiling_lib}_annotated{cutoff_suffix}_complete/bean_element_result.MultiMixtureNormal.csv",
        mixnormal_acc_res="results/model_runs/bean_negctrl/bean_run_result.bean_count_{tiling_lib}_annotated{cutoff_suffix}_complete/bean_element_result.MultiMixtureNormal+Acc.csv",
        crisphiermix_res_all="results/model_runs/CRISPhieRmix_negctrl/CRISPhieRmix_run_result.bean_count_{tiling_lib}_annotated_complete.target_allpyEdited/CRISPhieRmix.csv",
        crisphiermix_res_behive="results/model_runs/CRISPhieRmix_negctrl/CRISPhieRmix_run_result.bean_count_{tiling_lib}_annotated_complete.target_behive/CRISPhieRmix.csv",
        cb2_out_all="results/model_runs/CB2/CB2_run_result.bean_count_{tiling_lib}_annotated_complete.target_allEdited/CB2_gene.csv",
        cb2_out_pred="results/model_runs/CB2/CB2_run_result.bean_count_{tiling_lib}_annotated_complete.target_behive/CB2_gene.csv",
    params:
        bean_prefix="results/model_runs/bean_negctrl/bean_count_{tiling_lib}_annotated{cutoff_suffix}_complete/",
        mageck_prefix="results/model_runs/mageck_negctrl/bean_count_{tiling_lib}_annotated_complete/",
    output:
        "results/model_runs/bean_count_{tiling_lib}_annotated{cutoff_suffix}_complete/all_scores_negctrl.csv"
    run:
        shell("python scripts/evaluate_model_runs/collect_scores_tiling.py bean_count_{wildcards.tiling_lib}_annotated{wildcards.cutoff_suffix}_complete --result-suffix _negctrl --noallele-screen-name bean_count_{wildcards.tiling_lib}_annotated_complete")


rule run_2reps_tiling_negctrl:
    input:
        input_h5ad="results/filtered_annotated/{tiling_lib}/bean_count_{tiling_lib}_annotated{cutoff_suffix}_complete.h5ad"
    output:
        "results/model_runs/bean_count_{tiling_lib}_annotated{cutoff_suffix}_complete_rep3_rep4/all_scores_negctrl.csv"
    run:
        shell("python scripts/run_models/run_models_on_2_reps_tiling.py {input.input_h5ad} --use-negctrl ")