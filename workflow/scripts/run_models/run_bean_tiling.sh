bdata_path=$1
pids=$()

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.Normal.csv ] ; then
    bean-run tiling $bdata_path -o results/model_runs/bean --allele-df-key sig_allele_counts_spacer_2_7_translated_prop0.05_0.2 --cuda 
    pids+=($!)
fi

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.Normal.csv ] ; then
    bean-run tiling $bdata_path --scale-by-acc --acc-bw-path resources/accessibility/ENCFF262URW.hg19.bw -o results/model_runs/bean --allele-df-key sig_allele_counts_spacer_2_7_translated_prop0.05_0.2 --cuda 
    pids+=($!)
fi

for pid in ${pids[*]}; do
    wait $pid
done