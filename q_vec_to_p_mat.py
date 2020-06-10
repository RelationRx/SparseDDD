import numpy


def qvec_to_pmat(q_vec=None, model=None):
    start_vec = q_vec[-model.num_basis + 1:]
    vec_form = q_vec[1: -model.num_basis]
    q_mat = model.Indc.astype('float64')
    q_mat[model.Indc] = vec_form
    p_mat = numpy.linalg.solve(model.MassMat, q_mat)
    p_mat[:, -1] = start_vec
    # p_mat = full(p_mat)
    return p_mat
