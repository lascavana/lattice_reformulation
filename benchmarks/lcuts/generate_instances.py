import math
import numpy as np

def generate_aardal(filename, rng, n, m, type):
    # m constrains, n original variables and m slack variables
    assert type in ['binary', 'integer', 'unbounded']

    cost = rng.randint(low=1, high=1000, size=n)
    A = rng.randint(low=1, high=1000, size=(m,n))


    if type == 'integer':
        upper_bounds = rng.randint(low=5, high=10, size=(n,))
    else:
        upper_bounds = np.ones(n, dtype=int)

    rhs = 0.5*A.dot(upper_bounds)
    rhs = [math.floor(e) for e in rhs]

    with open(filename, 'w') as lp_file:
        lp_file.write("maximize\nOBJ:" + "".join([f" + {cost[j]} x{j+1}" for j in range(n)]) + "\n")
        lp_file.write("\nsubject to\n")

        for i in range(m):
            lp_file.write(f"C{i+1}:" + "".join([f" +{A[i,j]} x{j+1}" for j in range(n)]) + f" <= {rhs[i]}\n")

        # bounds of original variables
        if type == 'integer':
            lp_file.write("\nbounds\n")
            for j in range(n):
                lp_file.write(f"0 <= x{j+1} <= {upper_bounds[j]}\n")

        if type == 'binary':
            lp_file.write("\nbinary\n" + " ".join([f"x{j+1}" for j in range(n)]) + "\n")
        else:
            lp_file.write("\ngenerals\n" + " ".join([f"x{j+1}" for j in range(n)]) + "\n")



rng = np.random.RandomState(seed=0)
types = ['binary', 'integer', 'unbounded']

for type in types:
    m=1
    for n in [10, 20, 50, 100]:
        filename = f'{type}_m{m}_n{n}.lp'
        generate_aardal(filename, rng, n, m, type)
    n = 50
    for m in [2, 3, 4]:
        filename = f'{type}_m{m}_n{n}.lp'
        generate_aardal(filename, rng, n, m, type)
