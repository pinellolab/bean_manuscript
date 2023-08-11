bdata_path=$1
#control_label=$2
pids=()

echo "Submitting CRISPhieRmix runs..."

mkdir -p results/model_runs/CRISPhieRmix_negctrl/CRISPhieRmix_run_result.$(basename "${bdata_path}" .h5ad).target_{behive,allEdited}

Rscript scripts/run_models/run_CRISPhieRmix_LDLRCDS.R \
-i $bdata_path \
-o results/model_runs/CRISPhieRmix_negctrl/CRISPhieRmix_run_result.$(basename "${bdata_path}" .h5ad) &
pids+=($!)

for pid in ${pids[*]}; do
    wait $pid
done