import csv
import os
from pathlib import Path

if __name__ == '__main__':
    MY_PATH = Path(__file__)
    HOME_PATH = Path.home()
    EXEC_PATH = (MY_PATH / "../../../cmake-build-release/experiments").resolve()
    UCR_PATH = HOME_PATH / "Univariate_ts"
    EECSV_PATH = (MY_PATH / "../../eeOutputFold0.csv").resolve()
    RESULTS_OUTDIR_PATH = (MY_PATH / "../results/json").resolve()

    cmdsfile = open("cmdsfile", "w")

    print(f"SCRIPT FILE = {MY_PATH}")

    # Test if the exec path exists
    if not (EXEC_PATH.exists() and EXEC_PATH.is_file() and os.access(EXEC_PATH, os.X_OK)):
        print(f"Error: {EXEC_PATH} not found or is not an executable file")
        exit(1)
    print(f"EXEC_PATH = {EXEC_PATH}")

    # Test UCR Path
    if not (UCR_PATH.exists() and UCR_PATH.is_dir()):
        print(f"Error: {UCR_PATH} not found or not a directory")
        exit(1)
    print(f"UCR_PATH = {UCR_PATH}")

    # Test EECSV
    if not (EECSV_PATH.exists() and EECSV_PATH.is_file() and os.access(EECSV_PATH, os.R_OK)):
        print(f"Error: {EECSV_PATH} not found or is not a readable file")
        exit(1)
    print(f"EECSV_PATH = {EECSV_PATH}")

    # Test where the commands will write the results
    if not (RESULTS_OUTDIR_PATH.exists() and RESULTS_OUTDIR_PATH.is_dir()):
        print(f"Error: commands output directory {RESULTS_OUTDIR_PATH} does not exist or is not a directory")
        exit(1)
    print(f"RESULTS_OUTDIR_PATH = {RESULTS_OUTDIR_PATH}")



    # Generate all commands for a given running mode (when available)
    # Note: we do not deal with derivative here
    def generate_cmd(record, mode, output):
        # Extract all components from the record (obtained from EE CSV)
        name, \
        cdtw_w, \
        wdtw_weight, \
        __ddtw_w, \
        __wddtw_weight, \
        lcss_w, lcss_epsilon, \
        msm_cost, \
        twe_nu, twe_lambda, \
        erp_g, erp_w = record

        # Test if the folder exists
        folder = UCR_PATH / name
        if not (folder.exists() and folder.is_dir()):
            print(f"Error: {folder} does not exist")
            exit(1)

        head = f"{EXEC_PATH} -ucr {UCR_PATH}"

        # Generate commands
        def do_cdtw(lb, w):
            cmd = head+f" -dataset {name} {mode} -out {RESULTS_OUTDIR_PATH}/{name}_cdtw_{lb}_{mode}.json -dist cdtw {w} {lb}"
            print(cmd, file=output)

        def do_dtw(lb):
            cmd = head+f" -dataset {name} {mode} -out {RESULTS_OUTDIR_PATH}/{name}_dtw_{lb}_{mode}.json -dist dtw {lb}"
            print(cmd, file=output)

        # --- CDTW, window ratio, no lb
        do_cdtw("lb-none", cdtw_w)
        do_dtw("lb-none")
        # --- CDTW, window ratio, keogh
        do_cdtw("lb-keogh", cdtw_w)
        do_dtw("lb-keogh")
        # --- CDTW, window ratio, keogh2
        do_cdtw("lb-keogh2", cdtw_w)
        do_dtw("lb-keogh2")
        # --- CDTW, window ratio, webb
        do_cdtw("lb-webb", cdtw_w)
        do_dtw("lb-webb")

    # # # # # # # # # # # # # # # # # #

    with open(EECSV_PATH, newline='') as csvfile:
        records = csv.reader(csvfile)
        header = next(records)  # Skip header
        print(header)
        for r in records:
            for mode in ["base", "base_ea", "pru", "pru_la", "eap", "eap_la"]:
                generate_cmd(r, mode, cmdsfile)
