import csv
import pyscipopt
import numpy as np


n_instances = 5
problem = 'cuww'
result_file = f'results/{problem}.csv'
instance_path = f'benchmarks/{problem}'

formulations = ['ahl', 'orig', 'ahl*', 'pataki']
fieldnames = ['instance', 'seed', 'formulation', 'nnodes', 'time', 'status', 'gap']


## CREATE SCIP MODEL ##
tol = 1e-8
scip_parameters = {'limits/time': 3600,
                   'limits/memory': 4000,
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
                name = f'instance_{i+1}.lp'
                if formulation == "ahl": 
                    name = f'ahl_instance_{i+1}.lp'
                elif formulation == "ahl*": 
                    name = f'regahl_instance_{i+1}.lp'
                elif formulation == "pataki": 
                    name = f'pat_instance_{i+1}.lp'
                results = {'instance': f'instance_{i+1}',
                           'formulation': formulation,
                           'seed': seed }


                m.readProblem(f'{instance_path}/{name}')
                m.setParam('randomization/permutationseed', seed)
                m.setParams(scip_parameters)

                m.optimize()

                print("... solved")
                results['nnodes'] = m.getNTotalNodes()
                results['time'] = m.getSolvingTime()
                results['status'] = m.getStatus()
                results['gap'] = m.getGap()
                print(results)

                if m.getStatus() == 'optimal':
                    if optimalsol is None:
                        optimalsol = m.getObjVal()
                    else:
                        assert( abs(optimalsol - m.getObjVal()) < tol)

                m.freeProb()

                writer.writerow(results)
                csvfile.flush()