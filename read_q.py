import glob
import argparse
import numpy as np
import pandas as pd

def get_range1(valarray):
    maxval = np.amax(np.abs(valarray))
    minval = np.amin(np.abs(valarray))
    R = np.log10(max(1.0,maxval)) - np.log10(max(1.0,minval))
    return R

def get_range2(valarray):
    maxval = np.amax(valarray)
    minval = np.amin(valarray)

    if maxval == 0.0: maxval = 1.0
    if minval == 0.0: minval = 1.0

    R = np.sign(maxval)*np.log10(np.abs(maxval)) - np.sign(minval)*np.log10(np.abs(minval))
    return R

def get_Q_stats(filename):
    stats = {'filename': filename[2:]}

    # read file #
    with open(filename, 'r') as f:  
        lines = f.readlines()

    # get Q matrix #
    Q = []
    for line in lines:
        l = line.split(' ')[:-1]
        row = [int(val) for val in l]
        Q.append(row)
    Q = np.array(Q)

    m, n = Q.shape
    nonzero_x, nonzero_y = np.nonzero(Q)

    # length of the last column #
    stats['last'] = np.linalg.norm(Q[:,-1])

    # relative length last-to-first #
    stats['last_rel'] = np.linalg.norm(Q[:,-1]) / np.linalg.norm(Q[:,0])

    # matrix density # 
    stats['d'] = nonzero_x.shape[0] / (m * n)

    # max row density #
    values, counts = np.unique(nonzero_x, return_counts=True)
    ind = np.argmax(counts)
    max_row = int(values[ind])
    nonzero = np.nonzero(Q[max_row])[0]
    stats['rowd'] = nonzero.shape[0] / n

    return stats


def get_A_stats(filename):
    stats = {'filename': filename[5:]}

    with open(filename, 'r') as f:  
        lines = f.readlines()

    Aext = []
    for line in lines:
        l = line.split(' ')[:-1]
        row = [int(val) for val in l]
        Aext.append(row)
    Aext = np.array(Aext)
    A = Aext[:,:-1]
    b = Aext[:,-1]
    
    # range of A #
    range_A = get_range1(A)
    stats['range_A'] = range_A

    # range of b #
    range_b = get_range1(b)
    stats['range_b'] = range_b

    return stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
    'directory'
    )
    args = parser.parse_args()

    directory = args.directory
    filenames = glob.glob(f'{directory}/*.txt')

    # Aext stats #
    fields = ['filename', 'range_A', 'range_b']
    data = {field: [] for field in fields}

    for filename in filenames:
        stats = get_A_stats(filename)
        for field in fields:
            data[field].append(stats[field])


    # # Q stats #
    # fields = ['filename', 'last', 'last_rel', 'd', 'rowd']
    # data = {field: [] for field in fields}

    # for filename in filenames:
    #     stats = get_Q_stats(filename)
    #     for field in fields:
    #         data[field].append(stats[field])

    df = pd.DataFrame(data=data)
    df.to_excel('miplib.xlsx', sheet_name='sheet1', float_format="%.2f", index=False)





    