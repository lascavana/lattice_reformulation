import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from scipy.stats.mstats import gmean

def gmean(x):
    x = np.array(x)
    if np.where(x==0.0)[0].size > 0:
        return 0.0
    else:
        a = np.log(x)
        return np.exp(a.mean())

def MIPLIB_table():
    data = pd.read_csv("results/MIPLIB_subset_all.csv")
    results = {}

    instances = list(set(data['instance']))
    instances.sort()
    seeds = set(data['seed'])
    formulations = set(data['formulation'])

    for instance in instances:
        i_data = data.loc[data['instance'] == instance]
        nnodes = {'nnodes_ahl':[], 'nnodes_orig': []}
        for seed in seeds:
            is_data = i_data.loc[i_data['seed'] == seed]
            for formulation in formulations:
                isf_data = is_data.loc[is_data['formulation'] == formulation]
                nnodes[f'nnodes_{formulation}'].append((isf_data['nnodes'].values[0]))

        results[instance] = {'nnodes_ahl_mean': np.mean(nnodes['nnodes_ahl']),
                            'nnodes_orig_mean': np.mean(nnodes['nnodes_orig']),
                            'nnodes_ahl_gmean': gmean(nnodes['nnodes_ahl']),
                            'nnodes_orig_gmean': gmean(nnodes['nnodes_orig']),
                            'nnodes_ahl_std': np.std(nnodes['nnodes_ahl']),
                            'nnodes_orig_std': np.std(nnodes['nnodes_orig']),
                            'nnodes_ahl': nnodes['nnodes_ahl'],
                            'nnodes_orig': nnodes['nnodes_orig']}
    
        print(instance, 
              "{:.2f}".format(results[instance]['nnodes_orig_gmean']), 
              " & ",
              "{:.2f}".format(results[instance]['nnodes_ahl_gmean']) )



def box_plot(ax, data, edge_color, fill_color):
    bp = ax.boxplot(data, patch_artist=True)

    for element in ['boxes', 'whiskers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color)

    for patch in bp['boxes']:
        patch.set(facecolor=fill_color)

    for patch in bp['fliers']:
        patch.set(markeredgecolor=edge_color)

    return bp

csv_file = "results/msplit_feasibility/msplit_feasibility.csv"
txt_file = "results/determinants_mknapsack.txt"

# read csv file #
results = {}
data = pd.read_csv(csv_file)
instances = set(data['instance'])
seeds = set(data['seed'])
formulations = set(data['formulation'])

# read results according to metric of choice #
metric = 'nnodes'
for instance in instances:
    i_data = data.loc[data['instance'] == instance]
    vals = {}
    for formulation in formulations:
        vals[f'vals_{formulation}'] = []
    for seed in seeds:
        is_data = i_data.loc[i_data['seed'] == seed]
        for formulation in formulations:
            isf_data = is_data.loc[is_data['formulation'] == formulation]
            vals[f'vals_{formulation}'].append((isf_data[metric].values[0]))

    results[instance] = {}
    for fo in formulations:
        results[instance][f'vals_{fo}_mean'] = np.mean(vals[f'vals_{fo}'])
        results[instance][f'vals_{fo}_gmean'] = gmean(vals[f'vals_{fo}'])
        results[instance][f'vals_{fo}_std'] = np.std(vals[f'vals_{fo}'])
        results[instance][f'vals_{fo}'] = vals[f'vals_{fo}']


# print results per instance #
if False:
    formulations = ['orig', 'ahl', 'ahl*', 'pataki']
    s = f"              "
    for fo in formulations:
        s += f"{fo}  "
    print(s)
    for instance in instances:
        s = f"{instance} "
        for fo in formulations:
            s += f"{results[instance][f'vals_{fo}_gmean']:.2f}  "
        print(s)
    assert 0


# re-order instances #
orig_gmean = [results[instance]['vals_orig_gmean'] for instance in instances]
ind = np.argsort(orig_gmean)
orig = [results[instance]['vals_orig'] for instance in instances]
ahl = [results[instance]['vals_ahl'] for instance in instances]
pataki = [results[instance]['vals_pataki'] for instance in instances]
# orig = [orig[j] for j in ind]
# ahl = [ahl[j] for j in ind]
# pataki = [pataki[j] for j in ind]


# Creating plot
fig, ax = plt.subplots()
bp1 = box_plot(ax, orig, 'red', (1, 0, 0, .3))
bp2 = box_plot(ax, ahl, 'blue', (0, 0, 1, .3))
bp3 = box_plot(ax, pataki, 'green', (0, 1, 0, .3))
ax.legend([bp1["boxes"][0], bp2["boxes"][0], bp3["boxes"][0]], ['Original', 'AHL', 'KP'])

# show plot
fontsize = 20
ax.set_title('Market split', fontsize=fontsize)
ax.set_xlabel('Instance', fontsize=fontsize)
ax.set_ylabel('Number of nodes', fontsize=fontsize)
ax.set_xticks([])
ax.set_yscale('log')
# ax.set_ylim([50,1e7])
# plt.show()
plt.savefig('marketsplit.pdf')

