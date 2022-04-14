import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.mstats import gmean

data = pd.read_csv("results/results.csv")
results = {}

instances = set(data['instance'])
seeds = set(data['seed'])
formulations = set(data['formulation'])

for instance in instances:
    i_data = data.loc[data['instance'] == instance]
    nnodes = {'nnodes_ahl':[], 'nnodes_cea': [], 'nnodes_orig': []}
    for seed in seeds:
        is_data = i_data.loc[i_data['seed'] == seed]
        for formulation in formulations:
            idf_data = is_data.loc[is_data['formulation'] == formulation]
            nnodes[f'nnodes_{formulation}'].append((idf_data['nnodes'].values[0]))

    results[instance] = {'nnodes_ahl_mean': gmean(nnodes['nnodes_ahl']),
                         'nnodes_cea_mean': gmean(nnodes['nnodes_cea']),
                         'nnodes_orig_mean': gmean(nnodes['nnodes_orig']),
                         'nnodes_ahl_std': np.std(nnodes['nnodes_ahl']),
                         'nnodes_cea_std': np.std(nnodes['nnodes_cea']),
                         'nnodes_orig_std': np.std(nnodes['nnodes_orig']),}

for instance in instances:
    print(instance, results[instance]['nnodes_orig_mean'])

orig = [results[instance]['nnodes_orig_mean'] for instance in instances]
ahl = [results[instance]['nnodes_ahl_mean'] for instance in instances]
cea = [results[instance]['nnodes_cea_mean'] for instance in instances]

ind = np.argsort(orig)
orig = [orig[j] for j in ind]
ahl = [ahl[j] for j in ind]
cea = [cea[j] for j in ind]

plt.plot(orig, label='Original', marker='o', linestyle='')
plt.plot(ahl, label='AHL', marker='o', linestyle='')
plt.plot(cea, label='CEA', marker='o', linestyle='')
#plt.xlim(0,30)
plt.xlabel('Instance')
plt.ylabel('Number of nodes')
plt.semilogy()
#plt.ylim(0,10000)
plt.legend()
plt.show()
