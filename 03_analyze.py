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


def box_plot(ax, data, edge_color, fill_color):
    bp = ax.boxplot(data, patch_artist=True)

    for element in ['boxes', 'whiskers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color)

    for patch in bp['boxes']:
        patch.set(facecolor=fill_color)

    for patch in bp['fliers']:
        patch.set(markeredgecolor=edge_color)

    return bp


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'problem',
        choices=['struct_s', 'struct_b', 'nostruct_s', 'nostruct_b', 'ms', 'gap', 'ca'],
        )
    parser.add_argument(
        'output',
        choices=['per_instance', 'aggregated', 'boxplot'],
        )
    args = parser.parse_args()
    problem = args.problem
    output = args.output
    metric = 'nnodes'

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
            results[instance][fo] = gmean(if_data[metric].values.tolist())



    # OPTION1: results per instance and formulation #
    if output == 'per_instance':
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
        formulations = ['orig', 'ahl', 'ahl_poor', 'pataki']
        for fo in formulations:
            r = [results[instance][fo] for instance in instances]
            r = gmean(r)
            print(fo, f"{r:.2f}")


    # OPTION3: boxplot #
    if output == 'boxplot':
        if problem not in ['ca', 'gap', 'ms']:
            raise NotImplementedError
        
        results = {}
        for instance in instances:
            results[instance] = {}
            i_data = data.loc[data['instance'] == instance]
            for fo in formulations:
                if_data = i_data.loc[i_data['formulation'] == fo]
                results[instance][fo] = if_data[metric].values.tolist()

        orig = [results[instance]['orig'] for instance in instances]
        ahl = [results[instance]['ahl'] for instance in instances]
        ahl_poor = [results[instance]['ahl_poor'] for instance in instances]
        pataki = [results[instance]['pataki'] for instance in instances]

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

