import math
import numpy
import scipy.linalg as linalg
from calculate_g_matrix import calculate_g_matrix


def single_time_point_error(pf_mat=None, coeff_vec=None, start_vec=None, t_gap=None, model=None):
    if model.norm_fac == 'Relative':
        normalisation_fac = (coeff_vec.T * model.MassMat) * coeff_vec
    else:
        normalisation_fac = 1
    mat_g = linalg.expm(t_gap * pf_mat)
    result_vec = (mat_g * start_vec) - coeff_vec
    error = ((result_vec.T * model.MassMat) * result_vec) / normalisation_fac
    if math.isnan(error):
        error = math.nan
        grad_mat = numpy.empty(pf_mat.shape[1], pf_mat.shape[1] + 1)
        grad_mat[:] = numpy.NaN
    else:
        tmp_middle_mat = ((model.MassMat * result_vec) * start_vec.T)
        g_tmp = calculate_g_matrix(pf_mat, tmp_middle_mat, t_gap, normalisation_fac)
        g_tmp_2 = (((2 * mat_g.T) * model.MassMat) * result_vec) / normalisation_fac
        grad_mat = numpy.concat([g_tmp, g_tmp_2])
    return error, grad_mat
