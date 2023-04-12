bdata_path=$1
pids=()

echo "Submitting CB2 runs..."

if [ ! -f results/model_runs/CB2/CB2_run_result.$(basename "${bdata_path}" .h5ad)/CB2.csv ] ; then
    Rscript scripts/run_models/run_CB2.R \
    -i $bdata_path \
    -o results/model_runs/CB2/CB2_run_result.$(basename "${bdata_path}" .h5ad) &
    pids+=($!)
fi

for pid in ${pids[*]}; do
    wait $pid
done