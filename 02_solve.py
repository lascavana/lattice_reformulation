import csv
import pyscipopt
import numpy as np


n_instances = 1
problem = 'mknapsack'
result_file = f'results/{problem}.csv'
instance_path = f'benchmarks/{problem}'

formulations = ['ahl', 'orig']
fieldnames = ['instance', 'seed', 'formulation', 'nnodes', 'time', 'status']


## CREATE SCIP MODEL ##
tol = 1e-6
scip_parameters = {'limits/time': 3600,
                   'timing/clocktype': 1,
                   'numerics/feastol': tol,
                   'display/verblevel': 0}

m = pyscipopt.Model()


## SOLVE INSTANCES ##
with open(result_file, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for i in range(n_instances):
        optimalsol = None
        for seed in range(4):
            print(f"~~ Instance {i+1}, seed {seed}")
            for formulation in formulations:
                print(f"    Formulation {formulation}")
                name = f'instance_{i+1}_100_6.lp'
                if formulation == "ahl": name = f'ahl_instance_{i+1}_100_6.lp'
                results = {'instance': f'instance_{i+1}',
                           'formulation': formulation,
                           'seed': seed }


                m.readProblem(f'{instance_path}/{name}')
                m.setParam('randomization/permutationseed', seed)
                m.setParams(scip_parameters)

                m.optimize()

                print("... solved")
                results['nnodes'] = m.getNNodes()
                results['time'] = m.getSolvingTime()
                results['status'] = m.getStatus()
                print(results)

                if m.getStatus() == 'optimal':
                    if optimalsol is None:
                        optimalsol = m.getObjVal()
                    else:
                        assert( abs(optimalsol - m.getObjVal()) < tol)

                m.freeProb()

                writer.writerow(results)
                csvfile.flush()
