import numpy as np

def generate_markshare(filename, m, L, U, rng):
    n = 10*(m-1)
    A = rng.randint(low=L, high=U-1, size=(m,n))
    b = np.floor(np.sum(A, axis=1) / 2)
    b = b.astype(int)

    with open(filename, 'w') as lp_file:
        lp_file.write("minimize\nOBJ:" + "".join([f" +1 x{j+n+1}" for j in range(2*m)]) + "\n")
        lp_file.write("\nsubject to\n")

        for i in range(m):
            lp_file.write(f"C{i+1}:" + "".join([f" +{A[i,j]} x{j+1}" for j in range(n)]) + f" +1 x{n+1+2*i} -1 x{n+2+2*i}" + f" = {b[i]}\n")

        lp_file.write("\nbinary\n" + " ".join([f"x{j+1}" for j in range(n)]) + "\n")


m = 4
l = 0
u = 99
rng = np.random.RandomState(seed=0)

for k in range(30):
    generate_markshare(f'instance_{k+1}.lp', m, l, u, rng)
