import os
import matplotlib.pyplot as plt
import pathlib as path
import numpy as np
import json
import datetime
import math

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
            except:
                print("Skipping f")
    #
    return results


# From the doc https://matplotlib.org/examples/api/barchart_demo.html
def autolabel(ax, rects, fontsize):
    for rect in rects:
        height = rect.get_height()
        txt = f"{height:2.2f}"
        ax.text(rect.get_x() + rect.get_width() / 2., height / 2, txt, ha='center', va='bottom',
                fontsize=fontsize,
                color=(0, 0, 0),
                backgroundcolor=(0.9, 0.9, 0.9))


# Bar figure
def bar_fig(list_tuple, xlabel, ylabel, title, output_path):
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 16
    ax.tick_params(axis='both', labelsize=fontsize)
    x, val = zip(*list_tuple)
    bar = ax.bar(x, val)
    ax.set_xlabel(xlabel, fontsize=fontsize, labelpad=15)
    ax.set_ylabel(ylabel, fontsize=fontsize, labelpad=10)
    ax.set_title(title, fontsize=fontsize)
    autolabel(ax, bar, fontsize)
    fig.tight_layout()
    fig.show()
    fig.savefig(output_path + ".pdf", bbox_inches='tight', pad_inches=0)
    fig.savefig(output_path + ".png", bbox_inches='tight', pad_inches=0)
    plt.close(fig)


# Analyse per running mode
def analysis_per_running_mode(json_results, outdir):
    # Iterate over the results, aggregate by distance, mode, lb
    results_cdtw = {}  # Dictionary that will hold the aggregated results for CDTW
    results_dtw = {}  # Dictionary that will hold the aggregated results for DTW
    check_cdtw = {}  # Test if the results of nb correct are the same for all run
    check_dtw = {}
    modes = set()  # Set of discovered modes
    lbs = set()  # Set of discovered lower bound
    for r in json_results:
        # Check if the dataset is a "same length"
        train = r['train']
        test = r['test']
        ref = train['minl']
        if ref != train['maxl'] or ref != test['minl'] or ref != test['maxl']:
            print(f"Skipping {r['dataset']}")
            print("  ", train)
            print("  ", test)
            continue

        # switch dictionary according to the distance
        results = results_cdtw
        check = check_cdtw
        if r['distance']['name'] == "dtw":
            results = results_dtw
            check = check_dtw

        dataset_name = r['dataset']
        nb_correct = r['nb_correct']
        if not dataset_name in check:
            check[dataset_name] = nb_correct
        elif check[dataset_name] != nb_correct:
            print("Error: not matching result found for ")
            print(r)
            print("Expected ", check[dataset_name])
            exit(255)

        # aggregate by mode
        mode_name = r['distance']['runmode']
        modes.add(mode_name)
        mode_dic = results.setdefault(mode_name, {})

        # aggregate by lower bound
        lb_name = r['distance']['lb']
        lbs.add(lb_name)
        lb_dic = mode_dic.setdefault(lb_name, {})

        # sum time
        lb_total_runtime = lb_dic.setdefault('total_time_ns', int(0))
        runtime = r['timing_ns']
        total = lb_total_runtime + runtime
        lb_dic['total_time_ns'] = total
        lb_dic['total_time'] = str(datetime.timedelta(microseconds=total / 1e3))

    # --- Fix order of bounds
    lbs_order = ['lb-none', 'lb-keogh', 'lb-keogh2']
    if len(lbs_order) != len(lbs):
        print("Sort your bounds!")
        print("  Sorted: ", lbs_order)
        print("  Found:  ", lbs)
        exit(1)

    # --- Fix order of modes -- ignore lrp
    modes_order = ['base', 'pruneddtw', 'pru', 'base_ea', 'eap']
    modes_labels = {
        'base': "Base",
        'base_ea': "EABase",
        'pruneddtw': "PrunedDTW",
        'pru': "Pruned",
        'eap': "EAPruned"
    }
    modes_order_labels = [modes_labels[x] for x in modes_order]


    def analysis(results, distname, title):
        # --- Overall plot lbs/mode
        data_to_plot = []
        data_to_plot_hours = []
        for lb in lbs_order:
            modes_sorted = [(mode, results[mode][lb]['total_time_ns']) for mode in modes_order]
            in_minutes = [(m, float(ns) / (1e9 * 60)) for (m, ns) in modes_sorted]
            in_hours = [(m, float(ns) / (1e9 * 60 * 60)) for (m, ns) in modes_sorted]
            data_to_plot.append((lb, in_minutes))
            data_to_plot_hours.append((lb, in_hours))

        N = len(modes_order)
        M = len(lbs_order)
        fig = plt.figure(figsize=(30, 10))
        ax = fig.add_subplot(1, 1, 1)

        ind = np.arange(N)  # the x locations for the groups
        width = 0.23  # the width of the bars
        shift = 0
        group_bars = []
        group_lb = []
        for (lb, data) in data_to_plot_hours:
            minutes = [time for (mode, time) in data]
            p = ax.bar(ind + width * shift, minutes, width, bottom=0)
            shift += 1
            group_bars.append(p)
            group_lb.append(lb)

        fontsize = 28
        fontsize_label = 24
        ax.set_xlabel("lower bounds / mode", fontsize=fontsize, labelpad=15)
        ax.set_ylabel("runtime in hours", fontsize=fontsize, labelpad=10)
        ax.set_title(title, fontsize=fontsize)

        all_bars = [item for sublist in group_bars for item in sublist]
        autolabel(ax, all_bars, fontsize_label)

        ax.tick_params(axis='both', labelsize=fontsize)
        ax.set_xticks(ind + (len(data_to_plot) - 1) * width / 2)
        ax.set_xticklabels(modes_order_labels)

        ylabel = "runtime in hours"

        ax.legend(tuple([bg[0] for bg in group_bars]), tuple(group_lb), fontsize=fontsize)
        ax.autoscale_view()
        fig.tight_layout()

        output_path = os.path.join(outdir, f"{distname}")
        fig.savefig(output_path + ".pdf", bbox_inches='tight', pad_inches=0)
        fig.savefig(output_path + ".png", bbox_inches='tight', pad_inches=0)

    # --- ---

    # CDTW
    distname = "CDTW"
    title = f"{distname} runtime in minutes per mode and lower bound\n85 UCR Datasets, window parameter from EE"
    analysis(results_cdtw, distname, title)

    # DTW
    distname = "DTW"
    title = f"{distname} runtime in minutes per mode and lower bound\n85 UCR Datasets"
    analysis(results_dtw, distname, title)


# Scatter plot comparing mode,lb vs mode,lb
def datasets_plot(datasets_list, title, modex_name, lbx_name, modey_name, lby_name, outdir, pfx):
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(1, 1, 1)  # type:plt.axes.Axes

    # Get the data and create the scatter plot
    xs = []
    ys = []
    modex_name = modex_name.lower()
    modey_name = modey_name.lower()
    for d in datasets_list:
        modex = int(d[1][modex_name][lbx_name] / 1e6)  # convert from nano to milliseconds
        modey = int(d[1][modey_name][lby_name] / 1e6)
        xs.append(modex)
        ys.append(modey)
    sca = ax.scatter(x=xs, y=ys, s=5)

    # Adapt the ticks
    max_x = max(xs)
    max_y = max(ys)
    maxv = max(max_x, max_y)
    ticks = np.linspace(0, maxv, 6, dtype=int)
    scaling =int(10 ** (math.floor(math.log10(maxv))))
    mult = int(ticks[1] / scaling) * scaling
    if mult == 0:
        mult = scaling
    top = int((math.floor(maxv / mult) + 2) * mult)
    ticks = np.arange(0, top, mult)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.spines['right'].set_position(('outward', 0))

    # Square the range and draw the /.
    # Must be done after adapting the ticks!
    xrange = ax.get_xlim()
    yrange = ax.get_ylim()
    srange0 = int(min(xrange[0], yrange[0])) - 1
    srange1 = int(max(xrange[1], yrange[1])) + 1
    srange = (srange0, srange1)
    ax.set_xlim(srange)
    ax.set_ylim(srange)
    # Do the /
    orig0 = min(0, srange0)
    ax.plot([orig0, srange1], [orig0, srange1], "k", lw=0.5)

    # Use upper case names
    xname = modex_name + " " + lbx_name
    yname = modey_name + " " + lby_name
    xnameu = xname.upper()
    ynameu = yname.upper()

    # Auto text
    adiff = np.array(xs) - np.array(ys)
    xs_slower = np.sum(adiff > 0)
    xs_equals_ys = np.sum(adiff == 0)
    ys_slower = len(xs) - xs_slower - xs_equals_ys
    #print(xs_slower, xs_equals_ys, ys_slower)
    more_x = xnameu + f" is slower here: {xs_slower}"
    more_y = ynameu + f" is slower here: {ys_slower}"
    equal_xy = f"Draws: {xs_equals_ys}"
    # Create the text objects
    bot = min(srange)
    top = max(srange)
    base_style = dict(
        multialignment="left",
        horizontalalignment="center",
        verticalalignment="center",
        bbox=dict(facecolor='white', alpha=0.4)
    )
    ax.text(0.5 * top, 0.01 * top, more_x, **{**base_style, 'verticalalignment': "bottom"})
    ax.text(0.5 * top, 0.95 * top, more_y, **{**base_style, 'verticalalignment': "top"})
    ax.text(0.5 * top, 0.5 * top, equal_xy, **{**base_style, 'verticalalignment': "center"})

    # Labeling
    #
    title = title + "\n" + ynameu + " vs " + xnameu + " (in milliseconds)"
    ax.set_title(title, pad=20)
    #
    x_label = "Runtime for " + xnameu + " (in milliseconds)"
    ax.set_xlabel(x_label)
    #
    y_label = "Runtime for " + ynameu + " (in milliseconds)"
    ax.set_ylabel(y_label)
    # Use sci notation
    ax.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
    #
    fig.show()

    #ax.autoscale_view()
    fig.tight_layout()

    output_path = os.path.join(outdir, f"{pfx}{modex_name}_{lbx_name}_vs_{modey_name}_{lby_name}")
    fig.savefig(output_path + ".pdf", bbox_inches='tight', pad_inches=0)
    fig.savefig(output_path + ".png", bbox_inches='tight', pad_inches=0)


# Scatter plot Analysis
def analysis_per_dataset(json_results, outdir):
    # Analyse the results per distance / dataset / mode / lb
    results = {}
    for r in json_results:
        # Check if the dataset is a "same length"
        train = r['train']
        test = r['test']
        reflength = train['minl']
        if reflength != train['maxl'] or reflength != test['minl'] or reflength != test['maxl']:
            print(f"Skipping {r['dataset']}")
            print("  ", train)
            print("  ", test)
            continue
        #
        distance_name = r['distance']['name']
        all_datasets_dic = results.setdefault(distance_name, {})
        #
        dataset_name = r['dataset']
        dataset_dic = all_datasets_dic.setdefault(dataset_name, {})
        # Add info in the dictionary
        dataset_dic['series_length'] = reflength
        dataset_dic['train_size'] = train['size']
        dataset_dic['test_size'] = train['size']
        maxvalue = dataset_dic.setdefault('max', -math.inf)
        minvalue = dataset_dic.setdefault('min', math.inf)
        #
        mode_name = r['distance']['runmode']
        runmode_dic = dataset_dic.setdefault(mode_name, {})
        #
        lb_name = r['distance']['lb']
        runtime = r['timing_ns']
        runmode_dic[lb_name] = runtime
        dataset_dic['max'] = max(maxvalue, runtime)
        dataset_dic['min'] = min(minvalue, runtime)

    # Sort per max runtime, descending
    results_sorted = {}
    for (dist, all_datasets_dic) in results.items():
        l = sorted(all_datasets_dic.items(), key=lambda k: k[1]['max'], reverse=True)
        results_sorted[dist] = l

    # result_dist is a list of pair (dataset_name, json)
    def do_plots(result_dist, distname):
        lname = distname.lower()
        uname = distname.upper()

        # --- EAP
        datasets_plot(result_dist,
                      f"{uname} runtime comparison per dataset",
                      "eap", "lb-none", "eap", "lb-keogh", outdir, f"{lname}_")
        datasets_plot(result_dist,
                      f"{uname} runtime comparison per dataset",
                      "eap", "lb-none", "eap", "lb-keogh2", outdir, f"{lname}_")

        #datasets_plot(result_dist,
        #              f"{uname} runtime comparison per dataset",
        #              "eap_la", "lb-none", "ea", "lb-none", outdir, f"{lname}_")

        #datasets_plot(result_dist,
        #              f"{uname} runtime comparison per dataset",
        #              "eap_la", "lb-none", "eap", "lb-none", outdir, f"{lname}_")

        #datasets_plot(result_dist,
        #              f"{uname} runtime comparison per dataset",
        #              "eap_la", "lb-none", "eap_la_lrp", "lb-none", outdir, f"{lname}_")

        ## --- EAP LA + LB WEBB vs EAP LA + LB other
        #datasets_plot(result_dist,
        #              f"{uname} runtime comparison per dataset",
        #              "eap_la", "lb-webb", "eap_la", "lb-none", outdir, f"{lname}_")

        #datasets_plot(result_dist,
        #              f"{uname} runtime comparison per dataset",
        #              "eap_la", "lb-webb", "eap_la", "lb-keogh", outdir, f"{lname}_")

        #datasets_plot(result_dist,
        #              f"{uname} runtime comparison per dataset",
        #              "eap_la", "lb-webb", "eap_la", "lb-keogh2", outdir, f"{lname}_")


        # --- Sorting per series length size
        #dsize = sorted(result_dist, key=lambda k: k[1]['series_length'])
        #dhlen = int(len(dsize) / 2)
        #dshort = dsize[:dhlen]
        #dlong = dsize[dhlen:]
        #datasets_plot(dshort,
        #              f"{uname} runtime comparison for the 50% shortest datasets",
        #              "eap", "lb-none", "eap", "lb-keogh2", outdir, f"{lname}_short_")

        #datasets_plot(dlong,
        #              f"{uname} runtime comparison for the 50% longest datasets",
        #              "eap", "lb-none", "eap", "lb-keogh2", outdir, f"{lname}_long_")

    do_plots(results['dtw'].items(), 'dtw')
    do_plots(results['cdtw'].items(), 'cdtw')



if __name__ == "__main__":
    # Create the output directory
    outdir = "figures"
    if not os.path.exists(outdir): os.makedirs(outdir)

    # Read all the json
    json_results = read_all_in("results/json")

    # Create graph per running mode
    analysis_per_running_mode(json_results, outdir)

    # Create scatter plot
    #analysis_per_dataset(json_results, outdir)
