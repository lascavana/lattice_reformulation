import csv
import glob
import pyscipopt
import numpy as np

problem = 'MIPLIB_subset'
instances = glob.glob(f'benchmarks/{problem}/ahl_*.lp')
result_file = f'results/{problem}.csv'

formulations = ['ahl', 'orig']
fieldnames = ['instance', 'seed', 'formulation', 'nnodes', 'time', 'status', 'gap']


## CREATE SCIP MODEL ##
tol = 1e-6
isolate = False
scip_parameters = {'limits/time': 3600,
                   'timing/clocktype': 1,
                   'numerics/feastol': tol,
                   'display/verblevel': 0}
if isolate:
    scip_parameters.update(
        {'separating/maxroundsroot': 0,
         'separating/maxrounds': 0,
         'conflict/enable': FALSE,
         'conflict/useprop': FALSE,
         'presolving/maxrestarts': 0,
         'reoptimization/strongbranchinginit': FALSE}
    )
m = pyscipopt.Model()


## SOLVE INSTANCES ##
with open(result_file, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for instance in instances:
        optimalsol = None
        name = instance.split('/')[-1].split('ahl_')[-1].split('.')[0]

        for seed in range(4):
            print(f"~~ Instance {instance}, seed {seed}")
            for formulation in formulations:
                print(f"    Formulation {formulation}")
                results = {'instance': name,
                           'formulation': formulation,
                           'seed': seed }
                
                if formulation == "ahl": 
                    m.readProblem(instance)
                else:
                    m.readProblem(instance.replace("ahl_", "").replace(".lp",".mps"))

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
