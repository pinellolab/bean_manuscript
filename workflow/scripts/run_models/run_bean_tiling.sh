bdata_path=$1
pids=$()

if [[ "$bdata_path" == *"CBE"* ]]; then
    trail_args=" --allele-df-key sig_allele_counts_spacer_3_8_C.T_translated_prop0.05_0.2 --control-guide-tag CBE_CONTROL"
else
    trail_args=" --allele-df-key sig_allele_counts_spacer_2_7_translated_prop0.05_0.1 --control-guide-tag ABE_CONTROL"
fi

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.Normal_allEdited.csv ] ; then
    bean-run variant $bdata_path --perfect-edit -o results/model_runs/bean --target-col target_allEdited --cuda $trail_args &
    pids+=($!)
fi

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.Normal_behive.csv ] ; then
    bean-run variant $bdata_path --perfect-edit -o results/model_runs/bean --target-col target_behive --cuda $trail_args &
    pids+=($!)
fi

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.MultiMixtureNormal.csv ] ; then
    bean-run tiling $bdata_path -o results/model_runs/bean --splice-site-path results/model_runs/bean/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.MultiMixtureNormal.csv  --cuda $trail_args &
    pids+=($!)
fi

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.MultiMixtureNormal+Acc.csv ] ; then
    bean-run tiling $bdata_path --scale-by-acc --acc-bw-path resources/accessibility/ENCFF262URW.hg19.bw -o results/model_runs/bean --splice-site-path results/model_runs/bean/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.MultiMixtureNormal.csv --cuda $trail_args &
    pids+=($!)
fi

for pid in ${pids[*]}; do
    wait $pid
done