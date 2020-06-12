import numpy
import math
from python.single_time_point_error import single_time_point_error


def logical(array):
    return array != 0


def q_matrix_comparisons(q_mat=None, model=None, global_var=None, *args, **kwargs):
    start_vec = q_mat[-model.num_basis, -1]
    vec_form = q_mat[1, -1 - model.num_basis]
    q_mat = model.Indc
    q_mat[model.Indc] = vec_form
    error_ts = 0
    pf_mat = numpy.linalg.solve(model.MassMat, q_mat)
    # Full time
    if model.ShowComparison == 'FullSnapshot':
        err_mat = numpy.zeros(global_var.time_pts, 1)
        grad_mat = numpy.zeros(q_mat.shape[0], q_mat.shape[0] + 1)
        for vv in range(1, global_var.time_pts):
            t_gap = (global_var.time_pts(vv) - global_var.time_pts(1))
            coeff_vec = model.Coeff_vecs[vv]
            error_tmp, grad_mat_tmp = single_time_point_error(pf_mat, coeff_vec, start_vec, t_gap, model)
            err_mat[vv] = error_tmp
            grad_mat = grad_mat + grad_mat_tmp
        error = sum(err_mat) / global_var.num_tpts
        grad_mat = grad_mat / global_var.num_tpts
        grad_mat[:, range(1, model.num_basis)] = numpy.linalg.solve(model.MassMat, grad_mat[:, 1:model.num_basis])
        grad_mat = grad_mat(logical(numpy.concatenate([model.Indc, numpy.ones(model.num_basis, 1)])))

    # Evaluate first and last point only!
    if model.ShowComparison == 'FirstLastPoint':
        err_mat = numpy.zeros(2, 1)
        grad_mat = numpy.zeros(q_mat.shape[0], q_mat.shape[0] + 1)
        for vv in model.pts_of_int.reshape(-1):
            t_gap = (global_var.time_pts(vv) - global_var.time_pts(1))
            coeff_vec = model.Coeff_vecs[vv]
            error_tmp, grad_mat_tmp = single_time_point_error(pf_mat, coeff_vec, start_vec, t_gap, model)
            err_mat[vv] = error_tmp
            grad_mat = grad_mat + grad_mat_tmp
        err_mat = err_mat[model.pts_of_int]
        error = sum(err_mat) / 2
        grad_mat = grad_mat / 2
        grad_mat[:, 1: model.num_basis] = numpy.linalg.solve(model.MassMat, grad_mat[:, 1:model.num_basis])
        grad_mat = grad_mat(logical(numpy.concatenate([model.Indc, numpy.ones(model.num_basis, 1)])))

    # Time Series
    if model.ShowComparison == 'FullTimeSeries':
        err_mat = numpy.zeros(numpy.concatenate([model.TimeSeries.shape[0] - 1, model.n_ts_samples]))
        grad_mat = numpy.zeros(q_mat.shape[0], q_mat.shape[0] + 1)
        tot_valid = 0
        for ww in range(1, model.n_ts_samples):
            err_vec = numpy.zeros(1, global_var.time_pts - 1)
            grad_mat_cache = numpy.zeros(q_mat.size(), q_mat.shape[0] + 1)
            for vv in range(2, global_var.time_pts):
                t_gap = (global_var.time_pts(vv) - global_var.time_pts(1))
                coeff_vec = model.TimeSeries[vv, ww]
                error_tmp, grad_mat_tmp = single_time_point_error(pf_mat, coeff_vec, model.TimeSeries[1, ww],
                                                                  t_gap, model)
                err_vec[vv - 1] = error_tmp
                grad_mat_cache = grad_mat_cache + grad_mat_tmp
            if math.isnan(err_vec) == 0:
                err_mat[:, ww] = err_vec
                grad_mat = grad_mat + grad_mat_cache
                tot_valid = tot_valid + 1
        err_mat[:, sum(err_mat) == 0] = []
        error = sum(sum(err_mat)) / ((global_var.num_tpts - 1) * tot_valid)
        grad_mat = grad_mat / ((global_var.num_tpts - 1) * tot_valid)
        grad_mat[:, 1: model.num_basis] = numpy.linalg.solve(model.MassMat, grad_mat[:, 1: model.num_basis])
        grad_mat[:, -1] = 0
        grad_mat = grad_mat(logical(numpy.concat([model.Indc, numpy.ones(model.num_basis, 1)])))

    # Comparison
    if model.ShowComparison == 'Composite':
        if model.Lambda != 0:
            err_mat_ts = numpy.zeros(numpy.concatenate([model.TimeSeries.shape[0] - 1, model.n_ts_samples]))
            grad_mat_ts = numpy.zeros(q_mat.shape[0], q_mat.shape[0] + 1)
            tot_valid = 0
            for ww in range(1, model.n_ts_samples):
                err_vec = numpy.zeros(1, global_var.time_pts.size() - 1)
                grad_mat_cache = numpy.zeros(q_mat.shape[0], q_mat.shape[0] + 1)
                for vv in range(2, global_var.time_pts):
                    t_gap = (global_var.time_pts(vv) - global_var.time_pts(1))
                    coeff_vec = model.TimeSeries[vv, ww]
                    error_tmp, grad_mat_tmp = single_time_point_error(pf_mat, coeff_vec, model.TimeSeries[1, ww], t_gap,
                                                                      model)
                    err_vec[vv - 1] = error_tmp
                    grad_mat_cache = grad_mat_cache + grad_mat_tmp
                if math.isnan(err_vec):
                    err_mat_ts[:, ww] = err_vec
                    grad_mat_ts = grad_mat_ts + grad_mat_cache
                    tot_valid = tot_valid + 1
            err_mat_ts[:, sum(err_mat_ts) == 0] = []
            error_ts = sum(sum(err_mat_ts)) / ((global_var.num_tpts - 1) * tot_valid)
            grad_mat_ts = grad_mat_ts / ((global_var.num_tpts - 1) * tot_valid)
            grad_mat_ts[:, 1: model.num_basis] = numpy.linalg.solve(model.MassMat, grad_mat_ts[:, 1:model.num_basis])
            grad_mat_ts[:, -1] = 0
        # First last Snapshot
        err_mat_fl = numpy.zeros(2, 1)
        grad_mat_fl = numpy.zeros(q_mat.shape[0], q_mat.shape[0] + 1)
        for vv in model.pts_of_int.reshape(-1):
            t_gap = (global_var.time_pts(vv) - global_var.time_pts(1))
            coeff_vec = model.Coeff_vecs[vv]
            error_tmp, grad_mat_tmp = single_time_point_error(pf_mat, coeff_vec, start_vec, t_gap, model)
            err_mat_fl[vv] = error_tmp
            grad_mat_fl = grad_mat_fl + grad_mat_tmp
        err_mat_fl = err_mat_fl[model.pts_of_int]
        error_fl = sum(err_mat_fl) / 2
        grad_mat_fl = grad_mat_fl / 2
        grad_mat_fl[:, 1: model.num_basis] = numpy.linalg.solve(model.MassMat, grad_mat_fl[:, 1: model.num_basis])
        if model.Lambda != 0:
            error = (model.Lambda * error_ts) + ((1 - model.Lambda) * error_fl)
            grad_mat = (model.Lambda * grad_mat_ts) + ((1 - model.Lambda) * grad_mat_fl)
        else:
            error = error_fl
            grad_mat = grad_mat_fl
        err_mat = []
        grad_mat = grad_mat(logical(numpy.concatenate([model.Indc, numpy.ones(model.num_basis, 1)])))

    grad_mat[abs(grad_mat) > 10 ** 5] = 0
    # time keeping
    # stp_time = copy(toc)
    # if stp_time > 5:
    #     disp(strcat('slow eval:', num2str(stp_time)))

    return error, grad_mat, err_mat
