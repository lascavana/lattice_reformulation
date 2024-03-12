import os
import csv
import glob
import argparse
import pyscipopt
import numpy as np

def get_name(instance, formulation):
    problem = instance.split('/')[1]
    name = instance.split('/')[2]

    if formulation == 'orig':
        return name
    elif formulation == 'ahl':
        return 'ahl_' + name[:-3] + 'lp'
    elif formulation == 'pataki': 
        return 'pat_' + name[:-3] + 'lp'


if __name__ == "__main__":

    os.makedirs('results', exist_ok=True)
    result_file = 'results/MIPLIB.csv'
    instance_path = 'benchmarks/MIPLIB'

    formulations = ['orig', 'ahl', 'pataki']
    instances = glob.glob(f'{instance_path}/*.mps')
    fieldnames = ['instance', 'seed', 'formulation', 'nnodes', 'time', 'status', 'gap']

    # create PySCIPOpt model #
    tol = 1e-8
    scip_parameters = {'limits/time': 3600,
                        'limits/memory': 4000,
                        'timing/clocktype': 1,
                        'numerics/feastol': tol,
                        'display/verblevel': 0,
                        'randomization/permutevars': True}
    m = pyscipopt.Model()


    # solve instances #
    with open(result_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i, instance in enumerate(instances):
            for seed in range(5):
                print(f"~~ Instance {i+1}, seed {seed}")
                for formulation in formulations:
                    print(f"    Formulation {formulation}")

                    results = {'instance': instance,
                                'formulation': formulation,
                                'seed': seed }

                    name = get_name(instance, formulation)
                    m.readProblem(f'{instance_path}/{name}')
                    m.setParam('randomization/permutationseed', seed)
                    m.setParam('randomization/randomseedshift', seed)
                    m.setParams(scip_parameters)

                    m.optimize()

                    print("... solved")
                    results['nnodes'] = m.getNTotalNodes()
                    results['time'] = m.getSolvingTime()
                    results['status'] = m.getStatus()
                    results['gap'] = m.getGap()
                    print(results)

                    m.freeProb()

                    writer.writerow(results)
                    csvfile.flush()