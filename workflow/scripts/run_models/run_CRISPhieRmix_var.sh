bdata_path=$1
pids=()

echo "Submitting CRISPhieRmix runs..."

if [ ! -f results/model_runs/CRISPhieRmix/CRISPhieRmix_run_result.$(basename "${bdata_path}" .h5ad)/CRISPhieRmix.csv ] ; then
    Rscript scripts/run_models/run_CRISPhieRmix.R \
    -i $bdata_path \
    -o results/model_runs/CRISPhieRmix/CRISPhieRmix_run_result.$(basename "${bdata_path}" .h5ad) &
    pids+=($!)
fi

for pid in ${pids[*]}; do
    wait $pid
done