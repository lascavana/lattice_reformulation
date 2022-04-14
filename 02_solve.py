import csv
import pyscipopt
import numpy as np


n_instances = 2
formulations = ['ahl', 'orig']
instance_path = 'benchmarks/marksplit'
fieldnames = ['instance', 'seed', 'formulation', 'nnodes', 'time', 'status']


## CREATE SCIP MODEL ##
scip_parameters = {'limits/time': 3600,
                   'timing/clocktype': 1,
                   'numerics/feastol': 1e-7,
                   'display/verblevel': 0}

m = pyscipopt.Model()


## SOLVE INSTANCES ##
result_file = f"results.csv"
with open(result_file, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for i in range(n_instances):
        optimalsol = None
        for seed in range(4):
            print(f"~~ Instance {i}, seed {seed}")
            for formulation in formulations:
                print(f"    Formulation {formulation}")
                name = f'instance_{i+1}.lp'
                if formulation == "ahl": name = f'ahl_instance_{i+1}.lp'
                results = {'instance': f'instance_{i}',
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

                if m.getStatus() == 'optimal':
                    if optimalsol is None:
                        optimalsol = m.getObjVal()
                    else:
                        assert( abs(optimalsol - m.getObjVal()) < 1e-7)

                m.freeProb()

                writer.writerow(results)
                csvfile.flush()
