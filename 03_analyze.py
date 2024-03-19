import glob
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




if __name__ == "__main__":
    metric = 'nnodes'

    instances = glob.glob('benchmarks/MIPLIB/*.mps')
    instance_path = 'benchmarks/MIPLIB'
    instances = [f'{instance_path}/neos859080.mps']
    results = read_from_csv('MIPLIB', metric)

    formulations = ['ahl', 'orig']
    s = f"              "
    for fo in formulations:
        s += f"{fo}  "
    print(s)
    for instance in instances:
        s = f"{instance} "
        for fo in formulations:
            s += f"{results[instance][fo]:.2f}  & & "
        print(s)


