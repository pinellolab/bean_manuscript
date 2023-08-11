rule run_simulation:
    run:
        shell("python run_multi_simulation.py 20221222_acc_var simulated_screen_fullsort2 -n 10 -sa -sv")