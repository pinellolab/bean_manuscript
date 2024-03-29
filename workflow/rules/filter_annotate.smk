rule run_qc_notebooks:
    input:
        input_h5ad='results/mapped/{lib}/bean_count_{lib}_combined.h5ad'
    output:
        out_h5ad='results/filtered_annotated/{lib}/bean_count_{lib}_masked.h5ad',
        html_report='reports/qc_report_{lib}.html',
    params:
        html_prefix="reports/qc_report_{lib}"
    run:
        shell("bean-qc {input.input_h5ad} -o {output.out_h5ad} -r {params.html_prefix}")

rule filter_annotate_cds_alleles:
    input:
        input_h5ad='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_masked.h5ad',
        plasmid_h5ad='results/mapped/LDLRCDS/bean_count_LDLRCDS_plasmid.h5ad'
    params:
        output_prefix='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_alleleFiltered_spacer0.1_0.3'
    output:
        output_h5ad='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_alleleFiltered_spacer0.1_0.3.h5ad',
        output_filter_stats='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_alleleFiltered_spacer0.1_0.3.filtered_allele_stats.pdf',
        output_filter_log='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_alleleFiltered_spacer0.1_0.3.filter_log.txt',
    run:
        shell("bean-filter {input.input_h5ad} -o {params.output_prefix} -p {input.plasmid_h5ad} -s 0 -e 19 -w -b -t -ap 0.1 -sp 0.3")

rule filter_annotate_cds_alleles_CBE:
    input:
        input_h5ad='results/filtered_annotated/LDLRCDS_CBE_{cas_enzyme}/bean_count_LDLRCDS_CBE_{cas_enzyme}_masked.h5ad',
        plasmid_h5ad='results/mapped/LDLRCDS/bean_count_LDLRCDS_plasmid.h5ad'
    params:
        output_prefix='results/filtered_annotated/LDLRCDS_CBE_{cas_enzyme}/bean_count_LDLRCDS_CBE_{cas_enzyme}_alleleFiltered{cutoff_suffix}'
    output:
        output_h5ad='results/filtered_annotated/LDLRCDS_CBE_{cas_enzyme}/bean_count_LDLRCDS_CBE_{cas_enzyme}_alleleFiltered{cutoff_suffix}.h5ad',
        output_filter_stats='results/filtered_annotated/LDLRCDS_CBE_{cas_enzyme}/bean_count_LDLRCDS_CBE_{cas_enzyme}_alleleFiltered.filtered_allele_stats{cutoff_suffix}.pdf',
        output_filter_log='results/filtered_annotated/LDLRCDS_CBE_{cas_enzyme}/bean_count_LDLRCDS_CBE_{cas_enzyme}_alleleFiltered{cutoff_suffix}.filter_log.txt',
    run:
        shell("bean-filter {input.input_h5ad} -o {params.output_prefix} -p {input.plasmid_h5ad} -s 3 -e 8 -w -b -t -ap 0.1 -sp 0.3 ")

rule annotate_var:
    input:
        input_h5ad='results/filtered_annotated/LDLvar/bean_count_LDLvar_masked.h5ad'
    output:
        output_h5ad='results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated.h5ad',
    run:
        shell("python scripts/run_models/add_quantiles.py {input.input_h5ad} {output.output_h5ad}")

rule subset_complete_reps:
    input:
        input_h5ad='results/filtered_annotated/{library}/bean_count_{library}_annotated{cutoff_suffix}.h5ad'
    output:
        output_h5ad='results/filtered_annotated/{library}/bean_count_{library}_annotated{cutoff_suffix}_complete.h5ad',
    run:
        shell("python scripts/run_models/subset_complete_reps.py {input.input_h5ad} {output.output_h5ad}")

rule get_targetable_splice_pos:
    input:
        ldlr_fa="resources/LDLR/exons.fa"
    output:
        splice_targets="resources/LDLR/LDLR_ABE_splice_targets.csv"
    run:
        shell("python scripts/filter_annotate/get_splice_sites.py {input.ldlr_fa} A {output.splice_targets}")

rule annotate_tiling:
    input:
        input_h5ad='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_alleleFiltered{cutoff_suffix}.h5ad',
        behive_pred='resources/gRNA_info/target_prediction/LDLR-ABE_BEHive_consequence.xlsx',
        behive_ctrl_pred='resources/gRNA_info/target_prediction/control_gRNA_BEHive_consequence.csv',
        splice_sites='resources/LDLR/LDLR_ABE_splice_targets.csv'
    wildcard_constraints:
        cutoff_suffix="[_[a-z]*[.\d]+_[.\d]+]?"
    output:
        output_h5ad='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_annotated{cutoff_suffix}.h5ad',
    run:
        shell("python scripts/filter_annotate/assign_guide_to_outcome.py both {input.input_h5ad} {output.output_h5ad} -s {input.splice_sites} -p {input.behive_pred} --write-bdata --control-guide-tag ABE_CONTROL")
        shell("python scripts/run_models/add_quantiles.py {output.output_h5ad} {output.output_h5ad}")

rule annotate_tiling_CBEs:
    input:
        input_h5ad='results/filtered_annotated/LDLRCDS_CBE_{cas_enzyme}/bean_count_LDLRCDS_CBE_{cas_enzyme}_alleleFiltered{cutoff_suffix}.h5ad',
        behive_pred='resources/gRNA_info/target_prediction/LDLR-CBE_BEHive_consequence.csv',
        splice_sites='resources/LDLR/LDLR_CBE_splice_targets.csv'
    output:
        output_h5ad='results/filtered_annotated/LDLRCDS_CBE_{cas_enzyme}/bean_count_LDLRCDS_CBE_{cas_enzyme}_annotated{cutoff_suffix}.h5ad',
    run:
        shell("python scripts/filter_annotate/assign_guide_to_outcome.py all {input.input_h5ad} {output.output_h5ad} -s {input.splice_sites} -p {input.behive_pred} --write-bdata --control-guide-tag CBE_CONTROL")
        shell("python scripts/run_models/add_quantiles.py {output.output_h5ad} {output.output_h5ad}")


rule editing_pattern_analysis:
    input:
        var_h5ad='results/filtered_annotated/LDLvar/bean_count_LDLvar_masked.h5ad',
        cds_h5ad='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_masked.h5ad',
        notebook="notebooks/Fig1/1b_position.ipynb"
    output:
        pam_fig='notebooks/Fig1b/1b_pam_varcds_combined.pdf',
        pos_by_pam_ctxt_fig = 'notebooks/Fig1b/1b_pos_pam_and_context.pdf',
        pos_eff_behive_fig="notebooks/Fig1b/1b_pos_eff_behive_LDLRCDS.pdf",
        pos_eff_behive_normed_fig="notebooks/Fig1b/1b_pos_eff_behive_LDLRCDS_normed.pdf"
    run:
        shell("jupyter nbconvert --execute {input.notebook}")
        