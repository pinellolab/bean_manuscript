bdata_path=$1
pids=$()

if [[ "$bdata_path" == *"CBE"* ]]; then
    trail_args=( --control-guide-tag CBE_CONTROL --negctrl-col-value "CBE control" )
else
    trail_args=( --control-guide-tag ABE_CONTROL --negctrl-col-value "ABE control" )
fi

if [[ "$bdata_path" == *"CBE"* ]]; then
    trail_args2=( --allele-df-key sig_allele_counts_spacer_3_8_C.T_translated_prop0.05_0.2 )
else
    trail_args2=( --allele-df-key sig_allele_counts_spacer_0_19_A.G_translated_prop0.1_0.3 )
fi

if [[ "$bdata_path" == *"_complete"* ]]; then
    trail_args3=( --ignore-bcmatch )
else
    trail_args3=(  )
fi

if [ ! -f results/model_runs/bean_negctrl/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.Normal_allEdited.csv ] ; then
    bean-run variant $bdata_path --perfect-edit -o results/model_runs/bean_negctrl/ --splice-site-path results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.Normal.csv --target-column target_allEdited --fit-negctrl --negctrl-col Region --result-suffix _allEdited --cuda "${trail_args[@]}" "${trail_args3[@]}" &
    pids+=($!)
fi

if [ ! -f results/model_runs/bean_negctrl/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.Normal_behive.csv ] ; then
    bean-run variant $bdata_path --perfect-edit -o results/model_runs/bean_negctrl/ --splice-site-path results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.Normal.csv --target-column target_behive --fit-negctrl --negctrl-col Region --result-suffix _behive --cuda "${trail_args[@]}" "${trail_args3[@]}" &
    pids+=($!)
fi

if [ ! -f results/model_runs/bean_negctrl/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.MultiMixtureNormal.csv ] ; then
    bean-run tiling $bdata_path -o results/model_runs/bean_negctrl --splice-site-path results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.MultiMixtureNormal.csv  --fit-negctrl --negctrl-col Region --cuda "${trail_args[@]}" "${trail_args2[@]}" "${trail_args3[@]}" &
    pids+=($!)
fi

if [ ! -f results/model_runs/bean_negctrl/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.MultiMixtureNormal+Acc.csv ] ; then
    bean-run tiling $bdata_path --scale-by-acc --acc-bw-path resources/accessibility/ENCFF262URW.hg19.bw -o results/model_runs/bean_negctrl --splice-site-path results/model_runs/bean_negctrl/bean_run_result.bean_count_LDLRCDS_annotated/bean_element_result.MultiMixtureNormal+Acc.csv --fit-negctrl --negctrl-col Region --cuda "${trail_args[@]}" "${trail_args2[@]}" "${trail_args3[@]}" &
    pids+=($!)
fi

for pid in ${pids[*]}; do
    wait $pid
done