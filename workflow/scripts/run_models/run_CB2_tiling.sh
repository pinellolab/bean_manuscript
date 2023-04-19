bdata_path=$1
pids=()

echo "Submitting CB2 runs..."

for target_col in target_allEdited target_behive; do
    if [ ! -f results/model_runs/CB2/CB2_run_result.$(basename "${bdata_path}" .h5ad).$target_col/CB2_with_bcmatch_gene.csv ] ; then
        mkdir -p results/model_runs/CB2/CB2_run_result.$(basename "${bdata_path}" .h5ad).$target_col/
        Rscript scripts/run_models/run_CB2.R \
        -i $bdata_path \
        -o results/model_runs/CB2/CB2_run_result.$(basename "${bdata_path}" .h5ad).$target_col/ \
        -t $target_col &
        pids+=($!)
    fi
done

for pid in ${pids[*]}; do
    wait $pid
done