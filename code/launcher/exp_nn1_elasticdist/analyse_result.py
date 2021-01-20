import os
import matplotlib.pyplot as plt
import pathlib as path
import numpy as np
import json
import datetime

import pprint
pp = pprint.PrettyPrinter()


# Get the paths of all elements from a directory
def get_file_paths(directory):
    return [f for f in path.Path(directory).iterdir() if f.is_file()]

# Read all json file from a directory, return a list of dictionary
def read_all_in(directory):
    #
    results = []
    #
    all_file_paths = get_file_paths(directory)
    print(f"Found {len(all_file_paths)} files...")
    for fp in all_file_paths:
        with open(fp) as f:
            try:
                data = json.load(f)
                results.append(data)
            except json.JSONDecodeError as err:
                print("Cannot read file ", fp)
                print(err)
    #
    return results


# Bar figure
def bar_fig(list_tuple, xlabel, ylabel, title, output_path):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 16
    ax.tick_params(axis='both', labelsize=fontsize)
    x, val = zip(*list_tuple)
    bar = ax.bar(x, val)
    ax.set_xlabel(xlabel, fontsize=fontsize, labelpad=15)
    ax.set_ylabel(ylabel, fontsize=fontsize, labelpad=10)
    ax.set_title(title, fontsize=fontsize)
    # From the doc https://matplotlib.org/examples/api/barchart_demo.html
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            txt = f"{height:2.2f}"
            ax.text(rect.get_x() + rect.get_width() / 2., height / 2, txt, ha='center', va='bottom',
                    fontsize=fontsize,
                    color=(0, 0, 0),
                    backgroundcolor=(0.9, 0.9, 0.9))

    autolabel(bar)
    fig.tight_layout()
    fig.show()
    fig.savefig(output_path+".pdf", bbox_inches='tight', pad_inches=0)
    fig.savefig(output_path+".png", bbox_inches='tight', pad_inches=0)
    plt.close(fig)

if __name__ == "__main__":
    # Create the output directory
    outdir="figures"

    # Read all the json
    json_results = read_all_in("results/json")

    # Iterate over the results, aggregate by dataset, mode, distance
    results = {}  # Dictionary that will hold the aggregated results
    modes = set() # Set of discovered modes
    for r in json_results:
        # Aggregate by dataset
        # Note : setdefault is like get, but init to second arg if key is absent
        dataset_name = r['dataset']
        dataset_dic = results.setdefault(dataset_name, {})
        # Aggregate by mode. Sum runtime accross distances
        mode_name = r['distance']['runmode']
        modes.add(mode_name)
        mode_dict = dataset_dic.setdefault(mode_name, {})
        # Sum runtime across distances
        mode_total_runtime = mode_dict.setdefault('total_time_ns', int(0))
        runtime = r['timing_ns']
        total = mode_total_runtime + runtime
        mode_dict['total_time_ns'] = total
        mode_dict['total_time'] = str(datetime.timedelta(microseconds=total/1e3))
        # Distance
        distance_name = r['distance']['name']
        distance_dict = mode_dict.setdefault('distance', {})
        distance_dict[distance_name] = [runtime, str(datetime.timedelta(microseconds=runtime/1e3))]

    print(f"Loaded results for {len(results)} datasets")
    print(f"Found modes: {modes}")
    if len(results) == 0:
        print("No data, exit.")
        exit(0)

    # --- Total runtime per mode, without LCSS and SQED
    # Init
    modes_timings = {}
    for m in modes: modes_timings[m]=int(0)
    # Aggregate
    for _, dataset_dic in results.items():
        for m in modes:
            for k,v in dataset_dic[m]["distance"].items():
                if k == "lcss" or k == "sqed": continue
                modes_timings[m] += v[0]
    # Print
    modes_sorted = list(modes_timings.items())
    modes_sorted = sorted(modes_sorted, key=lambda x: x[1], reverse=True)
    modes_order  = list(zip(*modes_sorted))[0]
    for (m, t) in modes_sorted:
        print(f"{m:<10}: {datetime.timedelta(microseconds=t/1e3)}")

    # --- Graph
    doGraph = True
    if doGraph:
        if(not os.path.exists(outdir)): os.mkdir(outdir)

        # --- --- --- All distances
        in_hours = map(lambda x: (x[0], x[1]/(1e9*3600)), modes_sorted)
        xlabel = 'algorithm families'
        ylabel = 'runtime in hours'
        title =  "Total runtime in hours per algorithm family"\
                 +"\n85 UCR Datasets, parameters from EE"\
                 +"\nDTW, CDTW, WDTW, ERP, MSM, TWE"
        output = os.path.join(outdir, "all")
        bar_fig(in_hours, xlabel, ylabel, title, output)

        # --- --- --- Per distances
        per_dist = {}
        for (dataset_name, dataset) in results.items():
            for mode, timings in dataset.items():
                distances = timings["distance"]
                for dist, measures in distances.items():
                    d = per_dist.setdefault(dist, {})
                    runtime = d.setdefault(mode, 0)
                    d[mode] = runtime + measures[0]
        #pp.pprint(per_dist)
        for distname, dic in per_dist.items():
            k_v = []
            for m in modes_order:
                if m in dic:
                    k_v.append((m, dic[m]))
            in_hours = map(lambda x: (x[0], x[1] / (1e9 * 3600)), k_v)
            xlabel = 'algorithm families'
            ylabel = 'runtime in hours'
            title = "Total runtime in hours per algorithm family" \
                    + "\n85 UCR Datasets, parameters from EE" \
                    + "\n"+distname.upper()
            output = os.path.join(outdir, distname)
            bar_fig(in_hours, xlabel, ylabel, title, output)

        # --- Print per dataset
        doPerDataset = True
        if doPerDataset:
            # Runtime per dataset per mode
            for dataset_name, dataset_dic in results.items():
                print(dataset_name)
                l= []
                for m in modes_order:
                    acc=0
                    for k,v in dataset_dic[m]["distance"].items():
                        if k == "lcss" or k == "sqed": continue
                        acc += v[0]
                    l.append((m, acc, dataset_dic[m]['total_time']))
                for (m, _, time) in l: print(f"  {m:<10}: {time}")
                kv = map(lambda x: (x[0], x[1]/(1e9*3600)), l)
                xlabel = 'algorithm families'
                ylabel = 'runtime in hours'
                title = "Total runtime in hours per algorithm family" \
                        + "\n" + dataset_name + ", parameters from EE" \
                        + "\nDTW, CDTW, WDTW, ERP, MSM, TWE"
                output = os.path.join(outdir, "dataset_"+dataset_name)
                bar_fig(kv, xlabel, ylabel, title, output)




