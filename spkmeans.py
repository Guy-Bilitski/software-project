import sys
import numpy as np
import sys
import mykmeanssp

class Env:
    """ Class for global variables used in the system """
    k = "k"
    goal = "goal"
    input_file = "input_file"

def main():
    try:
        args = load_args()
    except Exception as ex:
        print(ex)
        return 1
    try:
        data_points = np.genfromtxt(fname=args.get(Env.input_file), dtype=float, delimiter=',')
        
        goal = args.get(Env.goal)
        if goal == 'spk':
            if args[Env.k] > data_points.shape[0]:
                print("Invalid input!")
                exit(1)
            T_as_pylist = mykmeanssp.transform_data_points(data_points.tolist(), args.get(Env.k))
            T_as_np = np.array(T_as_pylist)
            k = T_as_np.shape[1]
            indices, initial_centroids = kmeans_pp(T_as_np, k)
            final_centroids = mykmeanssp.kmeans(T_as_pylist, initial_centroids)
            print_centroids(indices, final_centroids)
        elif goal == 'wam':
            W = mykmeanssp.wam(data_points.tolist())
            print_matrix(W)
        elif goal == 'ddg':
            D = mykmeanssp.ddg(data_points.tolist())
            print_matrix(D)
        elif goal == 'lnorm':
            L = mykmeanssp.lnorm(data_points.tolist())
            print_matrix(L)
        elif goal == 'jacobi': # Sage - I removed lnorm here as the ouput correct using the input withput performing lnorm
            eigen_vectors, eigen_values = mykmeanssp.jacobi(data_points.tolist())
            print_jacobi_output(eigen_vectors, eigen_values)
        else:
            print("Invalid input!")
        
        #kmeans_pp(data_points, args.get(Env.k))
    except Exception as ex:
        print("An Error Has Occurred")
        return 1
    return 0



def kmeans_pp(data_points, k):
    try:
        indices, initial_centroids = get_centriods(data_points, k)
        return (indices, initial_centroids)
    except Exception as ex:
        raise Exception("An Error Has Occurred")

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
            
            
def print_matrix(matrix):
        for row in matrix:
            line = []
            for value in row:
                line.append('%.4f,' % value)
            if len(line) > 0:
                line[-1] = line[-1][:-1]
                print("".join(line))

def print_jacobi_output(matrix, eigenvalues):
    line = []
    for value in eigenvalues:
        line.append('%.4f,' % value)
    if len(line) > 0:
        line[-1] = line[-1][:-1]
        print("".join(line))
    
    for row in matrix:
        line = []
        for value in row:
            line.append('%.4f,' % value)
        if len(line) > 0:
            line[-1] = line[-1][:-1]
            print("".join(line))

def get_centriods(np_array, k):
    np.random.seed(0)
    n = np_array.shape[0]
    indices = [np.random.choice(n)]
    centroids = [np_array[indices[0], :]] # initializing the first centroid
    weighted_p = np.zeros(n, dtype=float)
    for _ in range(k - 1):
        for j in range(n):
            x = np_array[j,:]
            weighted_p[j] = (min(np.linalg.norm(x - c) for c in centroids))**2
        distance_sum = sum(weighted_p)
        np.divide(weighted_p, distance_sum, out=weighted_p)
        new_cent_index = np.random.choice(n, p=weighted_p)
        centroids.append(np.squeeze(np_array[new_cent_index, :]))
        indices.append(new_cent_index)
    centroids = [c.tolist() for c in centroids]
    indices = [int(i) for i in indices]
    return (indices, centroids)

def load_args():
    """ returns a dict with all args that mentioned in ENV class and inputted """
    inp = sys.argv

    args = {}
    if len(inp) != 4:
        raise Exception("Invalid Input!")
    try:
        args[Env.k] = int(inp[1])
        if args[Env.k] < 0:
            raise

        args[Env.goal] = inp[2]
        args[Env.input_file] = inp[3]
        
    except Exception:
        raise Exception("Invalid Input!")
    
    return args

main()