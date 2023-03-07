rule run_qc_notebooks:
    input:
        input_h5ad='results/mapped/{lib}/bean_count_{lib}_combined.h5ad'
    output:
        out_h5ad='results/filtered_annotated/{lib}/bean_count_{lib}_masked.h5ad',
        html_report='reports/qc_report_{lib}.html',
    params:
        html_prefix="reports/qc_report_{lib}"
    run:
        #shell("papermill notebooks/sample_quality_report.ipynb -p bdata_path {input.input_h5ad} -p out_bdata_path {output.out_h5ad} {output.report}")
        shell("bean-qc {input.input_h5ad} -o {output.out_h5ad} -r {params.html_prefix}")

rule filter_annotate_cds_alleles:
    input:
        input_h5ad='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_masked.h5ad',
        plasmid_h5ad='results/mapped/LDLRCDS/bean_count_plasmid_LDLRCDS.h5ad'
    params:
        output_prefix='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_alleleFiltered'
    output:
        output_h5ad='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_alleleFiltered.h5ad',
        output_filter_stats='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_alleleFiltered.filtered_allele_stats.pdf',
        output_filter_log='results/filtered_annotated/LDLRCDS/bean_count_LDLRCDS_alleleFiltered.filter_log.txt',
    run:
        shell("bean-filter {input.input_h5ad} {params.output_prefix} -p {input.plasmid_h5ad} -s 2 -s 7 -w -b -t -ap 0.05 ")

rule make_ldlvar_path:
    input:
        input_h5ad='results/filtered_annotated/LDLvar/bean_count_LDLvar_masked.h5ad'
    output:
        output_h5ad='results/filtered_annotated/LDLvar/bean_count_LDLvar_annotated.h5ad'
    run:
        shell("ln -s {input.input_h5ad} {output.output_h5ad}")
        