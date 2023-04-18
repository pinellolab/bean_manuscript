bdata_path=$1
pids=$()

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.MultiMixtureNormal.csv ] ; then
    bean-run tiling $bdata_path -o results/model_runs/bean --allele-df-key sig_allele_counts_spacer_2_7_translated_prop0.05_0.1 --splice-site-path results/model_runs/bean/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.MultiMixtureNormal.csv --control-guide-tag ABE_CONTROL 
    pids+=($!)
fi

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.MultiMixtureNormal+Acc.csv ] ; then
    bean-run tiling $bdata_path --scale-by-acc --acc-bw-path resources/accessibility/ENCFF262URW.hg19.bw -o results/model_runs/bean --allele-df-key sig_allele_counts_spacer_2_7_translated_prop0.05_0.1 --splice-site-path results/model_runs/bean/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.MultiMixtureNormal.csv --control-guide-tag ABE_CONTROL 
    pids+=($!)
fi

for pid in ${pids[*]}; do
    wait $pid
done