import numpy
from python.p_matrix_comparisons import p_matrix_comparisons
from python.q_matrix_comparisons import q_matrix_comparisons


def process_matrix(model, matrix, global_var, func):
    model.pts_of_int = (numpy.concatenate([1, 6]))
    a, b, _ = func(matrix, model, global_var, nargout=3)
    matrix[1] = matrix[1] + 0.0001
    A, B, _ = func(matrix, model, global_var, nargout=3)
    print(f'Derivative checking: {numpy.concatenate([(A - a) / 0.0001, b[1], B[1]])}')
    print(f'Error: {a}')
    return matrix


def process_p_matrix(model, p_mat, global_var):
    return process_matrix(model, p_mat, global_var, p_matrix_comparisons)


def process_q_matrix(model, q_mat, global_var):
    return process_matrix(model, q_mat, global_var, q_matrix_comparisons)
