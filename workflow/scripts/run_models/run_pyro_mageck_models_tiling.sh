prefix=$1
mask_col=$2
bdata_path=$3
mageck_result_prefix=$4
sub_reps=$5

rerun=false
pids=()
export OMP_NUM_THREADS=1
if [ $sub_reps != . ]; then
    sub_bdata_path=$(dirname "$bdata_path")/$( basename "${bdata_path}" .h5ad).${sub_reps//,/_}.h5ad
    bdata_path=$sub_bdata_path
    if [ ! -f $sub_bdata_path ];then
        echo "Subsetting replicates" $sub_reps "..."
        python subset_screen.py $bdata_path $sub_reps $sub_bdata_path
    fi
    echo "Using subsetted file at" $bdata_path
fi

if [ $mageck_result_prefix == . ]; then
    mageck_result_prefix=/PHShome/jr1025/projects/ANBE/mageck_results/
fi

echo "Submitting Pyro runs..."
echo "Trailing arguments" "${@:6}" "are passed to the python script."

for model in C16 C18 C19;#B18 B18_2;
#for model in A B16 B16_2 B16_0 B B2 B14 B14_2 B18 B18_2 B1;
do
    if [ ! -f ${prefix}/$(basename "${bdata_path}" .h5ad).model${model}.result.pkl ] ; then
        echo "python run_tiling_model_both.py -p=$prefix ${@:6} $bdata_path $model"
        if [ $mask_col != "." ] ; then
            python run_tiling_model_both.py -p=$prefix ${@:6} --sample-mask-column=$mask_col $bdata_path $model &
        else
            python run_tiling_model_both.py -p=$prefix ${@:6} $bdata_path $model &
        fi
    fi
    pids+=($!)
done
# MAGeCK input file /PHShome/jr1025/projects/ANBE/mageck_results/
for target_col in target_allEdited target_behive; do
    mageck_path=${mageck_result_prefix}/${prefix}/$(basename "${bdata_path}" .h5ad).$target_col/
    mkdir -p $mageck_path
    mageck_outfile=$mageck_path/$(basename "${bdata_path}" .h5ad).mageck_input.txt
    sgrna_eff=$mageck_path/$(basename "${bdata_path}" .h5ad).mageck_sgrna_eff.txt
    dm_topbot=${mageck_path}/$(basename "${bdata_path}" .h5ad).mageck_dm_topbot.txt
    dm_sort=${mageck_path}/$(basename "${bdata_path}" .h5ad).mageck_dm_sort.txt


    if [ $rerun ] || [ ! -f $mageck_outfile ] || [ ! -f $dm_topbot ] || [ ! -f $dm_sort ] || [ ! -f $sgrna_eff ]; then
        echo "Writing MAGeCK input files to " $mageck_path " ..."
        #if [[ $prefix == *"bcmatch"* ]] || [[ $prefix == *"lb"* ]] || [[ $prefix == *"sb"* ]]; then
        if [[ "$*" == *"use_bcmatch"* ]]; then
            echo "use barcode"
            mageck_outfile=$mageck_path/$(basename "${bdata_path}" .h5ad).bcmatch.mageck_input.txt
            if [ $mask_col != "." ] ; then
                python /PHShome/jr1025/projects/ANBE/mageck_results/make_mageck_input.py $bdata_path -p $mageck_path -m $mask_col --use_bcmatch --target_col $target_col
            else
                python /PHShome/jr1025/projects/ANBE/mageck_results/make_mageck_input.py $bdata_path -p $mageck_path  --use_bcmatch --target_col $target_col
            fi
        else
            if [ $mask_col != "." ] ; then
                python /PHShome/jr1025/projects/ANBE/mageck_results/make_mageck_input.py $bdata_path -p $mageck_path -m $mask_col --target_col $target_col
            else
                python /PHShome/jr1025/projects/ANBE/mageck_results/make_mageck_input.py $bdata_path -p $mageck_path --target_col $target_col 
            fi
        fi
    fi


    ## MAGeCK MLE
    echo "Submitting MAGeCK runs..."
    if [ $rerun ] || [ ! -f $mageck_path/topbot.gene_summary.txt ]; then
        mageck mle -k $mageck_outfile -d $dm_topbot -n $mageck_path/topbot --threads=1 &
        pids+=($!)
    fi
    if [ $rerun ] || [ ! -f $mageck_path/sort.gene_summary.txt ]; then
        mageck mle -k $mageck_outfile -d $dm_sort -n $mageck_path/sort --threads=1 &
        pids+=($!)
    fi
    if [ $rerun ] || [ ! -f $mageck_path/topbot_var.gene_summary.txt ]; then
        mageck mle -k $mageck_outfile -d $dm_topbot --genes-varmodeling 1000 -n $mageck_path/topbot_var --threads=1 &
        pids+=($!)
    fi
    if [ $rerun ] || [ ! -f $mageck_path/sort_var.gene_summary.txt ]; then
        mageck mle -k $mageck_outfile -d $dm_sort --genes-varmodeling 1000 -n $mageck_path/sort_var --threads=1 &
        pids+=($!)
    fi
done

for pid in ${pids[*]}; do
    wait $pid
done

# ## MAGeCK RRA
# bot_samples=$(awk -F'\t' '{if ($3=='0') {print $1}}' $dm_topbot | paste -sd "," -)
# top_samples=$(awk -F'\t' '{if ($3=='1') {print $1}}' $dm_topbot | paste -sd "," -)
# if [ $rerun ] || [ ! -f "$mageck_path/rra_top.gene_summary.txt" ]; then
#     echo "Running MAGeCK RRA..."
#     mageck test -k $mageck_outfile -t $top_samples -c $bot_samples --paired -n $mageck_path/rra_top 
# fi
# if [ $rerun ] || [ ! -f "$mageck_path/rra_bot.gene_summary.txt" ]; then
#     echo "Running MAGeCK RRA..."
#     mageck test -k $mageck_outfile -t $bot_samples -c $top_samples --paired -n $mageck_path/rra_bot 
# fi

# ## CRISPhieRmix

# ## ACE

# ## BAGEL