bdata_path=$1
pids=()

echo "Submitting CRISPRBean runs..."

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.Normal.csv ] ; then
    bean-run variant $bdata_path --uniform-edit -o results/model_runs/bean/ --load-existing &
    pids+=($!)
fi

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.MixtureNormal.csv ] ; then
bean-run variant $bdata_path -o results/model_runs/bean/ --load-existing &
pids+=($!)
fi

if [ ! -f results/model_runs/bean/bean_run_result.$(basename "${bdata_path}" .h5ad)/bean_element_result.MixtureNormal+Acc.csv ] ; then
bean-run variant $bdata_path --scale-by-acc --acc-bw-path resources/accessibility/ENCFF262URW.hg19.bw -o results/model_runs/bean --load-existing  &
pids+=($!)
fi

for pid in ${pids[*]}; do
    wait $pid
done