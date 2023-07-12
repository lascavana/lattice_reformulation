import os
import csv
import argparse
import pyscipopt
import numpy as np

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    'problem',
    choices=['struct_s', 'struct_b', 'nostruct_s', 'nostruct_b', 'ms', 'gap', 'ca'],
  )
  args = parser.parse_args()
  problem = args.problem
  n_instances = 30
  
  os.makedirs('results', exist_ok=True)
  result_file = f'results/{problem}.csv'
  instance_path = f'benchmarks/{problem}'

  formulations = ['orig', 'ahl', 'ahl_poor', 'pataki']
  if problem not in ['ms', 'gap', 'ca']:
    formulations.append('ahl_diag')

  fieldnames = ['instance', 'seed', 'formulation', 'nnodes', 'time', 'status', 'gap']


  ## create PySCIPOpt model ##
  tol = 1e-8
  scip_parameters = {'limits/time': 3600,
                    'limits/memory': 4000,
                    'timing/clocktype': 1,
                    'numerics/feastol': tol,
                    'display/verblevel': 0,
                    'randomization/permutevars': True}
  m = pyscipopt.Model()


  ## solve instances ##
  with open(result_file, 'w', newline='') as csvfile:
      writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
      writer.writeheader()

      for i in range(n_instances):
          optimalsol = None
          for seed in range(5):
              print(f"~~ Instance {i+1}, seed {seed}")
              for formulation in formulations:
                  print(f"    Formulation {formulation}")
                  name = f'instance_{i+1}.lp'
                  if formulation == "ahl": 
                      name = 'ahl_' + name
                  elif formulation == "ahl_diag": 
                      name = 'ahl_diag_' + name
                  elif formulation == "ahl_poor": 
                      name = 'ahl_poor_' + name
                  elif formulation == "pataki": 
                      name = 'pat_' + name
                  results = {'instance': f'instance_{i+1}',
                            'formulation': formulation,
                            'seed': seed }


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