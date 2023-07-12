import numpy as np

rng = np.random.RandomState(72)

structure = False
num_instances = 50
size = 100

for j in range(num_instances):
    
    if structure:
        filename = f'instance_{j+1}_n{size}_struct.txt'
        while True:
            M = rng.randint(1e4, 20001)
            N = rng.randint(1e3, 2001)
            p = rng.randint(1, 11, size)
            r = rng.randint(-10, 11, size)

            a = M*p + N*r

            if np.amin(a)>0: break
    else:
        filename = f'instance_{j+1}_n{size}_nostruct.txt'
        a = rng.randint(1e4, 150001, size)
    a = np.sort(a)

    with open(filename, 'w') as f:
        for num in a:
            f.write(str(num)+ " ")
