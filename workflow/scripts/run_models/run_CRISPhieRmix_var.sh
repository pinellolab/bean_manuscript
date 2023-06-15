bdata_path=$1
pids=()

if [ ! -f results/model_runs/CRISPhieRmix/CRISPhieRmix_run_result.$(basename "${bdata_path}" .h5ad)/CRISPhieRmix.csv ] ; then
    echo "Submitting CRISPhieRmix runs..."
    Rscript scripts/run_models/run_CRISPhieRmix.R \
    -i $bdata_path \
    -o results/model_runs/CRISPhieRmix/CRISPhieRmix_run_result.$(basename "${bdata_path}" .h5ad) &
    pids+=($!)
fi

for pid in ${pids[*]}; do
    wait $pid
done