set -e
bdata_path=$1

mask_col=mask
pids=()

mageck_result_prefix=results/model_runs/mageck


# MAGeCK input file /PHShome/jr1025/projects/ANBE/mageck_results/
for target_col in target_allEdited target_behive; do
    mageck_path=${mageck_result_prefix}/${prefix}/$(basename "${bdata_path}" .h5ad).$target_col/
    mkdir -p $mageck_path
    mageck_outfile=$mageck_path/$(basename "${bdata_path}" .h5ad).mageck_input.txt
    sgrna_eff=$mageck_path/$(basename "${bdata_path}" .h5ad).mageck_sgrna_eff.txt
    dm_topbot=${mageck_path}/$(basename "${bdata_path}" .h5ad).mageck_dm_topbot.txt
    dm_sort=${mageck_path}/$(basename "${bdata_path}" .h5ad).mageck_dm_sort.txt


    if [ ! -f $mageck_outfile ] || [ ! -f $dm_topbot ] || [ ! -f $dm_sort ] || [ ! -f $sgrna_eff ]; then
        echo "Writing MAGeCK input files to " $mageck_path " ..."
        mageck_outfile=$mageck_path/$(basename "${bdata_path}" .h5ad).bcmatch.mageck_input.txt
        python scripts/run_models/make_mageck_input.py $bdata_path -p $mageck_path -m $mask_col --use_bcmatch --target_col $target_col
    fi


    ## MAGeCK MLE
    echo "Submitting MAGeCK runs..."
    if [ ! -f $mageck_path/topbot.gene_summary.txt ]; then
        mageck mle -k $mageck_outfile -d $dm_topbot -n $mageck_path/topbot --threads=1 &
        pids+=($!)
    fi
    if [ ! -f $mageck_path/sort.gene_summary.txt ]; then
        mageck mle -k $mageck_outfile -d $dm_sort -n $mageck_path/sort --threads=1 &
        pids+=($!)
    fi
    if [ ! -f $mageck_path/topbot_var.gene_summary.txt ]; then
        mageck mle -k $mageck_outfile -d $dm_topbot --genes-varmodeling 1000 -n $mageck_path/topbot_var --threads=1 &
        pids+=($!)
    fi
    if [ ! -f $mageck_path/sort_var.gene_summary.txt ]; then
        mageck mle -k $mageck_outfile -d $dm_sort --genes-varmodeling 1000 -n $mageck_path/sort_var --threads=1 &
        pids+=($!)
    fi
done

for pid in ${pids[*]}; do
    wait $pid
done