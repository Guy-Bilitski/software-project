import sys
import pandas as pd
import numpy as np
import time
import mykmeanssp

class Env:
    """ Class for global variables used in the system """
    k = "k"
    epsilon = "epsilon"
    input_file1 = "input_file1"
    input_file2 = "input_file2"
    maxiter = "maxiter"

def main():
    try:
        args = load_args()
    except Exception as ex:
        print(ex)
        return 1
    try:
        input_data_frame = get_df(args.get(Env.input_file1),
                                  args.get(Env.input_file2))
        data_points = pd.DataFrame.to_numpy(input_data_frame, dtype=float)
        indices, initial_centroids = get_centriods(data_points, args.get(Env.k))
        data_points = [c.tolist() for c in data_points[:,1:]]
        final_centroids = mykmeanssp.fit(data_points, initial_centroids, args.get(Env.maxiter, -1), args.get(Env.epsilon))
        print_centroids(indices, final_centroids)
    except Exception as ex:
        print("An Error Has Occurred")
        return 1
    return 0


def print_centroids(indices, centroids):
    line = []
    for i in indices:
        line.append(f'{i},')
    if len(line) > 0:
        line[-1] = line[-1][:-1]
        print("".join(line))
    
    for c in centroids:
        line = []
        for i in c:
            line.append('%.4f,' % i)
        if len(line) > 0:
            line[-1] = line[-1][:-1]
            print("".join(line))


def get_centriods(np_array, k):
    np.random.seed(0)
    n = np_array.shape[0]
    indices = [np.random.choice(np_array[:,0])]
    centroids = [np.squeeze(np_array[:,1:][np_array[:,0]== indices[0]])] # initializing the first centroid
    weighted_p = np.zeros(n, dtype=float)
    for _ in range(k - 1):
        for j in range(n):
            x = np_array[j,1:]
            weighted_p[j] = (min(np.linalg.norm(x - c) for c in centroids))**2
        distance_sum = sum(weighted_p)
        np.divide(weighted_p, distance_sum, out=weighted_p)
        new_cent_index = np.random.choice(np_array[:,0], p=weighted_p)
        centroids.append(np.squeeze(np_array[:,1:][np_array[:,0] == new_cent_index]))
        indices.append(new_cent_index)
    centroids = [c.tolist() for c in centroids]
    indices = [int(i) for i in indices]
    return (indices, centroids)


def get_df(input_file1, input_file2):
    df1 = pd.read_csv(input_file1, header=None, dtype=float)
    df2 = pd.read_csv(input_file2, header=None, dtype=float)
    final_df = pd.merge(df1, df2, how='inner', on=0, copy=False, sort=True)
    return final_df


def load_args():
    """ returns a dict with all args that mentioned in ENV class and inputted """
    inp = sys.argv

    args = {}
    if len(inp) < 5 or len(inp) > 6:
        raise Exception("Invalid Input!")
    try:
        args[Env.k] = int(inp[1])
        if args[Env.k] < 1:
            raise

        args[Env.input_file1] = inp[-2]
        args[Env.input_file2] = inp[-1]
        epsilon = float(inp[-3])
        if epsilon <= 0:
            raise
        else:
            args[Env.epsilon] = epsilon

        if len(inp) == 6: # Checking maxiter exists
            maxiter = int(inp[2])
            if maxiter <= 0: # Validating maxiter value
                raise
            else:
                args[Env.maxiter] = maxiter # If exists, added to args.
        
    except Exception:
        raise Exception("Invalid Input!")
    
    return args

main()