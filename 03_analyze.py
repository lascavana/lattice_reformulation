import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def gmean(x):
    x = np.array(x)
    if np.where(x==0.0)[0].size > 0:
        return 0.0
    else:
        a = np.log(x)
        return np.exp(a.mean())


def read_from_csv(problem, metric, seedmean=True):
    # read csv file #
    csv_file = f'results/{problem}.csv'
    data = pd.read_csv(csv_file)
    instances = set(data['instance'])
    formulations = set(data['formulation'])

    # read results according to metric of choice #
    results = {}
    for instance in instances:
        results[instance] = {}
        i_data = data.loc[data['instance'] == instance]
        for fo in formulations:
            if_data = i_data.loc[i_data['formulation'] == fo]
            if seedmean:
                results[instance][fo] = gmean(if_data[metric].values.tolist())
            else:
                results[instance][fo] = if_data[metric].values.tolist()
    return results


def box_plot(ax, data, edge_color, fill_color):
    bp = ax.boxplot(data, patch_artist=True)

    for element in ['boxes', 'whiskers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color)

    for patch in bp['boxes']:
        patch.set(facecolor=fill_color)

    for patch in bp['fliers']:
        patch.set(markeredgecolor=edge_color)

    return bp


def get_reformulation_time(problem, formulation):
    if problem not in ['ca', 'gap', 'ms']:
        raise NotImplementedError

    with open(f'logs/{problem}_{formulation}.txt') as f:
        lines = f.readlines()

    times = []
    for line in lines:
        if line[:8] == "    took":
            line = line.split(' ')
            time = float(line[5])
            times.append(time)

    print(f'  The average reformulation time was: {np.mean(times):.0f} ms')
    print(f'  The maximum reformulation time was: {np.amax(times):.0f} ms')


def get_volumes(problem):
    if problem not in ['struct_s', 'struct_b', 'nostruct_s', 'nostruct_b']:
        raise NotImplementedError

    with open(f'logs/{problem}_ahl.txt') as f:
        lines = f.readlines()

    # original box #
    logvolumes = []
    for line in lines:
        if line[:12] == "    original":
            line = line.split(' ')
            logvol = float(line[7])
            logvolumes.append(logvol)
    print(f'  The average log-volume of the original box was: {np.mean(logvolumes):.2f}')

    # theoretical bound #
    logvolumes_1 = []
    logvolumes_2 = []
    for line in lines:
        if line[:15] == "    theoretical":
            line = line.split(' ')
            logvolumes_1.append(float(line[7]))
            logvolumes_2.append(float(line[9]))
    print(f'  The average analytical log-volume was: {np.mean(logvolumes_1):.2f} + {np.mean(logvolumes_2):.2f}')

    # ahl box #
    logvolumes = []
    for line in lines:
        if line[:17] == "    reformulation":
            line = line.split(' ')
            logvol = float(line[7])
            logvolumes.append(logvol)
    print(f'  The average log-volume of the AHL box was: {np.mean(logvolumes):.2f}')

    # ahl poor box #
    with open(f'logs/{problem}_ahl_poor.txt') as f:
        lines = f.readlines()
    logvolumes = []
    for line in lines:
        if line[:17] == "    reformulation":
            line = line.split(' ')
            logvol = float(line[7])
            logvolumes.append(logvol)
    print(f'  The average log-volume of the AHL_low box was: {np.mean(logvolumes):.2f}')

    # ahl diag box #
    with open(f'logs/{problem}_ahl_diag.txt') as f:
        lines = f.readlines()
    logvolumes = []
    for line in lines:
        if line[:17] == "    reformulation":
            line = line.split(' ')
            logvol = float(line[7])
            logvolumes.append(logvol)
    print(f'  The average log-volume of the AHL_diag box was: {np.mean(logvolumes):.2f}')

    # kp box #
    with open(f'logs/{problem}_kp.txt') as f:
        lines = f.readlines()
    logvolumes = []
    for line in lines:
        if line[:17] == "    reformulation":
            line = line.split(' ')
            logvol = float(line[7])
            logvolumes.append(logvol)
    print(f'  The average log-volume of the KP box was: {np.mean(logvolumes):.2f}')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'problem',
        choices=['struct_s', 'struct_b', 'nostruct_s', 'nostruct_b', 'ms', 'gap', 'ca'],
        )
    parser.add_argument(
        'output',
        choices=['per_instance', 'aggregated', 'boxplot', 'reftime'],
        )
    args = parser.parse_args()
    problem = args.problem
    output = args.output
    metric = 'nnodes'


    # OPTION1: results per instance and formulation #
    if output == 'per_instance':
        results = read_from_csv(problem, metric)

        formulations = ['orig', 'ahl', 'ahl_poor', 'pataki']
        s = f"              "
        for fo in formulations:
            s += f"{fo}  "
        print(s)
        for j in range(30):
            instance = f'instance_{j+1}'
            s = f"{instance} "
            for fo in formulations:
                s += f"{results[instance][fo]:.2f}  & & "
            print(s)


    # OPTION2: aggregated results per formulation #
    if output == 'aggregated':
        results = read_from_csv(problem, metric)

        formulations = ['orig', 'ahl', 'ahl_poor', 'pataki']
        for fo in formulations:
            r = [results[instance][fo] for instance in results.keys()]
            r = gmean(r)
            print(fo, f"{r:.2f}")


    # OPTION3: boxplot #
    if output == 'boxplot':
        if problem not in ['ca', 'gap', 'ms']:
            raise NotImplementedError
        
        results = read_from_csv(problem, metric, seedmean=False)

        orig = [results[instance]['orig'] for instance in results.keys()]
        ahl = [results[instance]['ahl'] for instance in results.keys()]
        ahl_poor = [results[instance]['ahl_poor'] for instance in results.keys()]
        pataki = [results[instance]['pataki'] for instance in results.keys()]

        # Creating plot
        fig, ax = plt.subplots()
        bp1 = box_plot(ax, orig, 'red', (1, 0, 0, .3))
        bp2 = box_plot(ax, ahl, 'blue', (0, 0, 1, .3))
        bp3 = box_plot(ax, pataki, 'green', (0, 1, 0, .3))
        bp4 = box_plot(ax, ahl_poor, 'orange', (0.9, .8, .4, .3))
        ax.legend([bp1["boxes"][0], bp2["boxes"][0], bp4["boxes"][0], bp3["boxes"][0]], ['Original', 'AHL', 'AHL_low', 'KP'])

        # show plot
        fontsize = 20
        titles = {'ca': 'Combinatorial Auctions',
                  'gap': 'Generalized Assignment',
                  'ms': 'Market Split'}
        ax.set_title(titles[problem], fontsize=fontsize)
        ax.set_xlabel('Instance', fontsize=fontsize)
        ax.set_ylabel('Number of nodes', fontsize=fontsize)
        ax.set_xticks([])
        ax.set_yscale('log')
        # ax.set_ylim([50,1e7])
        # plt.show()
        plt.savefig(f'results/{problem}.pdf')


    # OPTION4: reformulation time #
    if output == 'reftime':
        for fo in ['ahl', 'ahl_poor', 'kp']:
            print(f'Formulation: {fo}')
            get_reformulation_time(problem, fo)

