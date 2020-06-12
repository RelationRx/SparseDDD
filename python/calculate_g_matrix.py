import numpy


def calculate_g_matrix(pf_mat=None, tmp_middle_mat=None, t_gap=None, normalisation_fac=None):
    p_ft = pf_mat.T
    tol = 1e-10
    s_k0 = tmp_middle_mat
    s_k1 = (p_ft * s_k0) + (s_k0 * p_ft)
    # Zeroth team
    scale_fac_tmp = t_gap / 1
    tmp_grad_mat = (s_k0 * scale_fac_tmp)
    # First term
    scale_fac_tmp = (t_gap * scale_fac_tmp) / 2
    tmp_grad_mat = tmp_grad_mat + (s_k1 * scale_fac_tmp)
    trigger = 0
    aa = 2
    # Second term onwards.
    while trigger == 0:

        s_k2 = (p_ft * s_k1) + (s_k1 * p_ft) - ((p_ft * s_k0) * p_ft)
        scale_fac_tmp = (t_gap * scale_fac_tmp) / (aa + 1)
        mat_to_add = (s_k2 * scale_fac_tmp)
        if numpy.any(numpy.isinf(mat_to_add)):
            break
        tmp_grad_mat = tmp_grad_mat + mat_to_add
        s_k0 = s_k1
        s_k1 = s_k2
        if numpy.max(mat_to_add) < (tol * normalisation_fac):
            trigger = aa
        else:
            aa = aa + 1

    return (2 * tmp_grad_mat) / normalisation_fac
