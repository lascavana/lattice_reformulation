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


MIPLIB_table()
assert 0
data = pd.read_csv("results/results_mknapsack_nontrans.csv")
results = {}

instances = set(data['instance'])
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



orig_mean = [results[instance]['nnodes_orig_mean'] for instance in instances]
ind = np.argsort(orig_mean)

orig = [results[instance]['nnodes_orig'] for instance in instances]
ahl = [results[instance]['nnodes_ahl'] for instance in instances]
determinants = [determinants[instance] for instance in instances]
orig = [orig[j] for j in ind]
ahl = [ahl[j] for j in ind]
determinants = [determinants[j] for j in ind]
print(determinants)

fig, ax = plt.subplots()

# Creating plot
bp1 = box_plot(ax, orig, 'red', (1, 0, 0, .3))
bp2 = box_plot(ax, ahl, 'blue', (0, 0, 1, .3))
ax.legend([bp1["boxes"][0], bp2["boxes"][0]], ['Original', 'Reformulated'])

# show plot
plt.title('Multiple knapsack')
plt.xticks([])
plt.xlabel('Instance')
plt.yscale('log')
plt.ylabel('Number of nodes')
plt.show()




## Same but with time ##
# instances = set(data['instance'])
# seeds = set(data['seed'])
# formulations = set(data['formulation'])
#
# for instance in instances:
#     i_data = data.loc[data['instance'] == instance]
#     time = {'time_ahl':[], 'time_orig': []}
#     for seed in seeds:
#         is_data = i_data.loc[i_data['seed'] == seed]
#         for formulation in formulations:
#             isf_data = is_data.loc[is_data['formulation'] == formulation]
#             time[f'time_{formulation}'].append((isf_data['time'].values[0]))
#
#     results[instance] = {'time_ahl_mean': np.mean(time['time_ahl']),
#                          'time_orig_mean': np.mean(time['time_orig']),
#                          'time_ahl_gmean': gmean(time['time_ahl']),
#                          'time_orig_gmean': gmean(time['time_orig']),
#                          'time_ahl_std': np.std(time['time_ahl']),
#                          'time_orig_std': np.std(time['time_orig']),
#                          'time_ahl': time['time_ahl'],
#                          'time_orig': time['time_orig']}
#
#
#
# orig_mean = [results[instance]['time_orig_mean'] for instance in instances]
# ind = np.argsort(orig_mean)
#
# orig = [results[instance]['time_orig'] for instance in instances]
# ahl = [results[instance]['time_ahl'] for instance in instances]
# orig = [orig[j] for j in ind]
# ahl = [ahl[j] for j in ind]
#
# fig, ax = plt.subplots()
#
# # Creating plot
# bp1 = box_plot(ax, orig, 'red', (1, 0, 0, .3))
# bp2 = box_plot(ax, ahl, 'blue', (0, 0, 1, .3))
# ax.legend([bp1["boxes"][0], bp2["boxes"][0]], ['Original', 'Reformulated'])
#
# # show plot
# plt.title('Multiple knapsack')
# plt.xticks([])
# plt.xlabel('Instance')
# plt.yscale('log')
# plt.ylabel('Time')
# plt.show()
