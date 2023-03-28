rule run_bean_var:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated.h5ad"
    output:
        normal_res="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated/bean_element_result.Normal.csv",
        mixnormal_res="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated/bean_element_result.MixtureNormal.csv",
        mixnormal_acc_res="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated/bean_element_result.MixtureNormal+Acc.csv",
    run:
        shell("sh scripts/run_models/run_bean_var.sh {input.input_h5ad}")

rule run_mageck_var:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated.h5ad"
    output:
        mle_out="results/model_runs/mageck/bean_count_LDLvar_annotated/sort.gene_summary.txt",
        mle_var_out="results/model_runs/mageck/bean_count_LDLvar_annotated/sort_var.gene_summary.txt",
        mle_em_out="results/model_runs/mageck/bean_count_LDLvar_annotated/EM/sort.gene_summary.txt",
        mle_var_em_out="results/model_runs/mageck/bean_count_LDLvar_annotated/EM/sort_var.gene_summary.txt",
        mle_emf_out="results/model_runs/mageck/bean_count_LDLvar_annotated/EMf/sort.gene_summary.txt",
        mle_var_emf_out="results/model_runs/mageck/bean_count_LDLvar_annotated/EMf/sort_var.gene_summary.txt",
        rra_out="results/model_runs/mageck/bean_count_LDLvar_annotated/rra_top.gene_summary.txt",
        rra_bot_out="results/model_runs/mageck/bean_count_LDLvar_annotated/rra_bot.gene_summary.txt",
    run:
        shell("sh scripts/run_models/run_mageck_var.sh {input.input_h5ad} .")

rule evaluate_varscreen:
    input:
        mle_out="results/model_runs/mageck/bean_count_LDLvar_annotated/sort.gene_summary.txt",
        mle_var_out="results/model_runs/mageck/bean_count_LDLvar_annotated/sort_var.gene_summary.txt",
        mle_em_out="results/model_runs/mageck/bean_count_LDLvar_annotated/EM/sort.gene_summary.txt",
        mle_var_em_out="results/model_runs/mageck/bean_count_LDLvar_annotated/EM/sort_var.gene_summary.txt",
        mle_emf_out="results/model_runs/mageck/bean_count_LDLvar_annotated/EMf/sort.gene_summary.txt",
        mle_var_emf_out="results/model_runs/mageck/bean_count_LDLvar_annotated/EMf/sort_var.gene_summary.txt",
        rra_out="results/model_runs/mageck/bean_count_LDLvar_annotated/rra_top.gene_summary.txt",
        rra_bot_out="results/model_runs/mageck/bean_count_LDLvar_annotated/rra_bot.gene_summary.txt",
        normal_res="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated/bean_element_result.Normal.csv",
        mixnormal_res="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated/bean_element_result.MixtureNormal.csv",
        mixnormal_acc_res="results/model_runs/bean/bean_run_result.bean_count_LDLvar_annotated/bean_element_result.MixtureNormal+Acc.csv",
    params:
        bean_prefix="results/model_runs/bean/bean_count_LDLvar_annotated",
        mageck_prefix="results/model_runs/mageck/bean_count_LDLvar_annotated/",
    output:
        "results/model_runs/bean_count_LDLvar_annotated/splicing.metrics.xlsx",
        "results/model_runs/bean_count_LDLvar_annotated/strongest_splicing.metrics.xlsx",
    run:
        shell("python scripts/evaluate_model_runs/get_performance_varscreen.py bean_count_LDLvar_annotated")
        shell("python scripts/evaluate_model_runs/get_performance_varscreen.py bean_count_LDLvar_annotated --control strongest_splicing")

rule run_2reps_varscreen:
    input:
        input_h5ad="results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated.h5ad"
    output:
        "results/model_runs/bean_count_LDLvar_annotated_rep14_rep15/all_scores.csv"
    run:
        shell("python scripts/run_models/run_models_on_2_reps.py {input.input_h5ad}")

rule run_bean_tiling_acc:
    input:
        input_h5ad="results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_annotated.h5ad",
        splice_site="resources/LDLR/LDLR_ABE_splice_targets.csv"
    output:
        res="results/model_runs/bean/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.MultiMixtureNormal+Acc.csv",
    run:
        shell("bean-run tiling {input.input_h5ad} --scale-by-acc --acc-bw-path resources/accessibility/ENCFF262URW.hg19.bw -o results/model_runs/bean --allele-df-key sig_allele_counts_spacer_2_7_translated_prop0.05_0.1 --splice-site-path {input.splice_site} --control-guide-tag ABE_CONTROL ")

rule run_bean_tiling:
    input:
        input_h5ad="results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_annotated.h5ad",
        splice_site="resources/LDLR/LDLR_ABE_splice_targets.csv"
    output:
        res="results/model_runs/bean/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.MultiMixtureNormal.csv",
    run:
        shell("bean-run tiling {input.input_h5ad} -o results/model_runs/bean --allele-df-key sig_allele_counts_spacer_2_7_translated_prop0.05_0.1 --splice-site-path {input.splice_site} --control-guide-tag ABE_CONTROL ")

rule run_mageck_tiling:
    input:
        input_h5ad="results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_annotated.h5ad"
    output:
        mle_out_all="results/model_runs/mageck/bean_count_LDLRCDS_annotated.target_allEdited/sort.gene_summary.txt",
        mle_out_pred="results/model_runs/mageck/bean_count_LDLRCDS_annotated.target_behive/sort.gene_summary.txt",
        mle_var_out_all="results/model_runs/mageck/bean_count_LDLRCDS_annotated.target_allEdited/sort_var.gene_summary.txt",
        mle_var_out_pred="results/model_runs/mageck/bean_count_LDLRCDS_annotated.target_behive/sort_var.gene_summary.txt",
    run:
        shell("sh scripts/run_models/run_mageck_tiling.sh {input.input_h5ad} .")

rule evaluate_tiling:
    input:
        mle_out_all="results/model_runs/mageck/bean_count_LDLRCDS_annotated.target_allEdited/sort.gene_summary.txt",
        mle_out_pred="results/model_runs/mageck/bean_count_LDLRCDS_annotated.target_behive/sort.gene_summary.txt",
        mle_var_out_all="results/model_runs/mageck/bean_count_LDLRCDS_annotated.target_allEdited/sort_var.gene_summary.txt",
        mle_var_out_pred="results/model_runs/mageck/bean_count_LDLRCDS_annotated.target_behive/sort_var.gene_summary.txt",
        # rra_out="results/model_runs/mageck/bean_count_LDLRCDS_annotated/rra_top.gene_summary.txt",
        # rra_bot_out="results/model_runs/mageck/bean_count_LDLRCDS_annotated/rra_bot.gene_summary.txt",
        mixnormal_res="results/model_runs/bean/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.MultiMixtureNormal.csv",
        mixnormal_acc_res="results/model_runs/bean/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.MultiMixtureNormal+Acc.csv",
    params:
        bean_prefix="results/model_runs/bean/bean_count_LDLRCDS_annotated",
        mageck_prefix="results/model_runs/mageck/bean_count_LDLRCDS_annotated/",
    output:
        "results/model_runs/bean_count_LDLRCDS_annotated/allEdited_negctrl.metrics.csv",
        "results/model_runs/bean_count_LDLRCDS_annotated/allEdited_syn.metrics.csv",
        "results/model_runs/bean_count_LDLRCDS_annotated/behive_negctrl.metrics.csv",
        "results/model_runs/bean_count_LDLRCDS_annotated/behive_syn.metrics.csv",
        "results/model_runs/bean_count_LDLRCDS_annotated/clinvar_pb.metrics.csv",
        "results/model_runs/bean_count_LDLRCDS_annotated/clinvar_plpb.metrics.csv",
    run:
        shell("python scripts/evaluate_model_runs/get_performance_tiling.py bean_count_LDLRCDS_annotated")

rule run_2reps_tiling:
    input:
        input_h5ad="results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_annotated.h5ad"
    output:
        "results/model_runs/bean_count_LDLRCDS_annotated_rep7_rep8/all_scores.csv"
    run:
        shell("python scripts/run_models/run_models_on_2_reps_tiling.py {input.input_h5ad} ")