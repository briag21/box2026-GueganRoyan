from python_tsp.heuristics import solve_tsp_local_search
from python_tsp.exact import solve_tsp_dynamic_programming
import numpy
from hw3_q1 import overlap


def init_distance_matrix(S):
    distance_mat = numpy.zeros((len(S) + 1, len(S) + 1))
    distance_mat[0][0] = numpy.inf
    for i in range(len(S)):
        distance_mat[0][i + 1] = len(S[i])
        distance_mat[i + 1][0] = 0

        for j in range(len(S)):
            if i != j:
                distance_mat[i + 1][j + 1] = len(S[j]) - overlap(S[i], S[j])
            else:
                distance_mat[i + 1][i + 1] = numpy.inf

    return distance_mat


def tspToSCS(S, exact):
    """if exact:=1 exact method, heuristic method if 0"""
    S = list(S)
    distance_mat = init_distance_matrix(S)
    if exact:
        h_path, _ = solve_tsp_dynamic_programming(
            distance_mat
        )  # the cost of the path is not needed in our case
    else:
        h_path, _ = solve_tsp_local_search(distance_mat)

    h_openpath = [i - 1 for i in h_path[1:]]

    s_star = S[h_openpath[0]]
    prev_i = h_openpath[0]

    for i in h_openpath[1:]:
        s_star += S[i][overlap(S[prev_i], S[i]) :]
        prev_i = i

    return s_star
