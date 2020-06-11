import math

import numpy
from numpy.linalg import norm
from scipy import sparse
from scipy.io import loadmat
from scipy.spatial.distance import cdist
from scipy.special import gamma
from sklearn.cluster import k_means
import ghalton
from mc_integration import mc_integration
from p_matrix_comparisons import p_matrix_comparisons
from q_matrix_comparisons import q_matrix_comparisons
from q_vec_to_p_mat import qvec_to_pmat
from utils import process_p_matrix, process_q_matrix

if __name__ == "__main__":
    # Load processed data
    Model = {}
    Global = loadmat('SDE_ProcessedData.mat')
    Global = Global['Global']
    # Rescale time to (0,1) interval
    Global["time_pts"] = Global["time_pts"] / (Global["time_pts"][-1] - Global["time_pts"][0])
    Model.norm_fac = 'Relative'
    # Choose number and locations of basis function
    Model.num_basis = 30
    _, Model.Basis_Loc = k_means(Global["all_pts_MixModel"], Model["num_basis"], 'Replicates', 5)
    Model.subfunc = {}
    # Define form of basis functions used.
    Model.subfunc.An = lambda n=None, eps=None: 2 * math.pi ** (n / 2) / gamma(n / 2)
    Model.subfunc.Vn = lambda n=None, eps=None: (eps ** n) * (math.pi ** (n / 2)) / gamma((n / 2) + 1)
    Model.subfunc.Intn = lambda n=None, eps=None: eps ** n / (n * (n + 1))
    Model.basis_nd = lambda x=None, n=None, eps=None: max(0, 1 - (norm(x.T) / eps)) / (
        (Model.subfunc.An(n, eps) * Model.subfunc.Intn(n, eps)))
    # Choose scale parameter for each basis function. Connect to (minimum) N nearest basis functions.
    Model.subfunc.connections = 2

    Model.Basis_eps = (numpy.zeros(Model.num_basis, 1))
    for ii in range(1, Model.num_basis):
        # print(ii)
        d_tmp = cdist(Model.Basis_Loc[ii, :], Model.Basis_Loc)
        d_tmp = numpy.sort(d_tmp[d_tmp > 0])
        Model.Basis_eps[ii] = 1.00001 * d_tmp[Model.subfunc.connections]

    # Use Halton Node Placing for Monte Carlo integration
    sequencer = ghalton.Halton(Global.Data_Dim)
    p = haltonset(Global.Data_Dim)
    X0 = net(p, 10 ** 5)
    for ii in range(1, Model.num_basis):
        # # print(ii)
        for jj in range(ii, Model.num_basis):
            tmp = mc_integration(ii, jj, Model, Global, X0)
            Model.MassMat[ii, jj] = tmp
            Model.MassMat[jj, ii] = tmp

    Model.MassMat = sparse.csr_matrix(Model.MassMat)
    # Plot Mass Matrix and Inverse
    Model.InvMassMat = numpy.linalg.inv(Model.MassMat)
    Model.InvMassMat = (Model.InvMassMat + Model.InvMassMat.T) / 2
    # Find coefficient vectors for each time pt.
    Model.Coeff_vecs = cell(Global.time_pts.size(), 1)
    for vv in range(1, Global.time_pts):
        Model.Coeff_vecs[vv] = numpy.zeros(Model.num_basis, 1)
        for ii in range(1, Model.num_basis).reshape(-1):
            tmp = Model.basis_nd(Global.sample_pts[vv] - Model.Basis_Loc[ii, :], Global.Data_Dim,
                                 Model.Basis_eps(ii))
            Model.Coeff_vecs[vv][ii] = sum(tmp)
        Model.Coeff_vecs[vv] = Model.Coeff_vecs[vv] / sum(Model.Coeff_vecs[vv])

    #
    Model.TimeSeries = cell(Global.time_pts.size(), Global.Num_pts_each_time[1])
    for ww in range(1, Global.Num_pts_each_time(1)):
        # print(ww)
        for vv in range(1, Global.time_pts):
            Model.TimeSeries[vv, ww] = numpy.zeros(Model.num_basis, 1)
            for ii in range(1, Model.num_basis):
                tmp = Model.basis_nd(Global.sample_pts[vv][ww, :] - Model.Basis_Loc[ii, :],
                                     Global.Data_Dim, Model.Basis_eps(ii))
                Model.TimeSeries[vv, ww][ii] = tmp
            Model.TimeSeries[vv, ww] = Model.TimeSeries[vv, ww] / sum(Model.TimeSeries[vv, ww])

    # Plot graphs of a few time points.
    # G_tmp = graph(Model.MassMat)
    # G_tmp = rmedge(G_tmp, range(1, numnodes(G_tmp)), range(1, numnodes(G_tmp)))

    #  Define useful quantities for later.
    Model.InvSums = sum(Model.InvMassMat)
    Model.Indc = sparse(Model.MassMat != 0)
    Indc_sum = sum(Model.Indc)
    Indc_zero = sparse(Model.MassMat == 0)
    Indc_vec = find(Model.MassMat != 0)
    Indc_diag = find(ismember(Indc_vec.T, range(1, Model.num_basis ** 2, (Model.num_basis + 1))))
    # Define initial guess at P matrix
    init_P = double(Model.Indc)
    for vv in range(1, Model.num_basis).reshape(-1):
        init_P[vv, vv] = 0
        init_P[vv, vv] = - sum(init_P[:, vv])

    # Calculate corresponding Q and vectorise non-zero entries
    init_Q = (Model.MassMat * init_P)
    init_Q[Indc_zero] = 0
    init_Q_vec = init_Q[logical_not(Indc_zero)]
    # Specify effective dimension of problem.
    eff_var = init_Q_vec.size()
    # # print(strcat('Effective number of parameters:', num2str(eff_var)))
    # # print(strcat('Eff. basis functions for full DDD:', num2str(sqrt(eff_var))))
    # Append initial condition onto end of vector.
    init_Q_vec[range(-1 + 1, -1 + Model.num_basis)] = Model.Coeff_vecs[1]
    # Creating constraints matrix for Aeq * init_Q_vec = beq
    Q_Aeq = sparse(Model.num_basis + 1, Model.num_basis + eff_var)
    run_tot_tmp = 1
    for kk in range(1, Model.num_basis).reshape(-1):
        ind_tmp = range(run_tot_tmp, run_tot_tmp + Indc_sum(kk) - 1)
        Q_Aeq[kk, ind_tmp] = Model.InvSums(Model.Indc[kk, :])
        run_tot_tmp = run_tot_tmp + Indc_sum(kk)

    # clear('kk', 'ind_tmp', 'run_tot_tmp')
    Q_Aeq[Model.num_basis + 1, range((-1 - Model.num_basis + 1), -1)] = 1
    Q_beq = sparse(Model.num_basis + 1, 1)
    Q_beq[-1] = 1
    # imagesc(numpy.concatenate(
    #     [numpy.concatenate([[sum(dot(Model.InvMassMat, init_Q)).T], [nan]]), dot(Q_Aeq, init_Q_vec), Q_beq]))
    # Create positivity constraint Alt * init_Q_vec <= blt
    Q_Alt = sparse(Model.num_basis + eff_var, Model.num_basis + eff_var)
    run_tot_tmp = 1
    for kk in range(1, Model.num_basis).reshape(-1):
        xx = find(Model.Indc[kk, :])
        yy = Model.InvMassMat(xx, xx)
        ind_tmp = range(run_tot_tmp, run_tot_tmp + Indc_sum(kk) - 1)
        Q_Alt[ind_tmp, ind_tmp] = - yy
        run_tot_tmp = run_tot_tmp + Indc_sum(kk)

    ind_tmp = range(run_tot_tmp, run_tot_tmp + Model.num_basis - 1)
    Q_Alt[ind_tmp, ind_tmp] = - numpy.eye(Model.num_basis)
    for kk in Indc_diag.reshape(-1):
        Q_Alt[kk, :] = - Q_Alt[kk, :]

    Q_blt = sparse(Model.num_basis + eff_var, 1)
    # clear('kk', 'ind_tmp')
    ##

    opts = optimoptions('fmincon', 'Display', 'iter-detailed')
    opts.UseParallel = true(1, 1)
    opts.MaxIterations = 3000
    opts.SpecifyObjectiveGradient = true(1, 1)
    opts.Algorithm = 'sqp'
    opts.MaxIterations = 300
    test_init_Q_vec = init_Q_vec
    test_init_Q_vec[range(1, eff_var)] = 0
    Model.ShowComparison = 'FullSnapshot'
    Model.Guess_Q_FS, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), test_init_Q_vec, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.Guess_Q_FS, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), Model.Guess_Q_FS, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.Guess_Q_FS, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), Model.Guess_Q_FS, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.ShowComparison = ('FullSnapshot')
    Outputs = {}
    Outputs.Q = {}
    # Model.Guess_Q_FS(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
    Outputs.Q.mean_err_FS, Grad_calc, Outputs.Q.err_mat_FS = q_matrix_comparisons(Model.Guess_Q_FS, Model,
                                                                                  Global, nargout=3)
    # print('=========================================================')
    # print('Best we can hope for:')
    # print(Outputs.Q.err_mat_FS)
    # print(strcat('mean error:', num2str(mean(Outputs.Q.err_mat_FS))))
    ##
    Model.ShowComparison = 'FirstLastPoint'
    Model.pts_of_int = numpy.concatenate([1, 4])
    Model.Guess_Q_FL, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), test_init_Q_vec, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.Guess_Q_FL, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), Model.Guess_Q_FL, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.Guess_Q_FL, _, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), Model.Guess_Q_FL, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.ShowComparison = ('FullSnapshot')
    # Model.Guess_Q_FL(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
    Outputs.Q.mean_err_FL, Grad_calc, Outputs.Q.err_mat_FL = q_matrix_comparisons(Model.Guess_Q_FL, Model,
                                                                                  Global, nargout=3)
    # print('=========================================================')
    # print('Worst case scenario:')
    # print(Outputs.Q.err_mat_FL)
    # print(strcat('mean error:', num2str(mean(Outputs.Q.err_mat_FL))))
    ##
    Model.ShowComparison = 'FullTimeSeries'
    Model.n_ts_samples = 100
    Model.Guess_Q_TS, _, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), test_init_Q_vec, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.Guess_Q_TS, _, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), Model.Guess_Q_TS, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.Guess_Q_TS, _, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), Model.Guess_Q_TS, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.ShowComparison = 'FullSnapshot'
    Model.Guess_Q_TS[range(-1 - Model.num_basis + 1, -1)] = Model.Coeff_vecs[1]
    Outputs.Q.mean_err_TS, Grad_calc, Outputs.Q.err_mat_TS = q_matrix_comparisons(Model.Guess_Q_TS, Model,
                                                                                  Global, nargout=3)
    # print('=========================================================')
    # print('Full Time Series:')
    # print(Outputs.Q.err_mat_TS)
    # print(strcat('mean error:', num2str(mean(Outputs.Q.err_mat_TS))))
    ##
    Model.ShowComparison = 'Composite'
    Model.Lambda = 0
    Model.n_ts_samples = 10
    Model.pts_of_int = numpy.concatenate([1, 4])
    Model.Guess_Q_C, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), test_init_Q_vec, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    ##
    Model.ShowComparison = 'Composite'
    Model.Lambda = 0.01
    Model.n_ts_samples = 10
    Model.pts_of_int = numpy.concatenate([1, 4])
    Model.Guess_Q_C, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), Model.Guess_Q_C, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    ##
    Model.ShowComparison = 'Composite'
    Model.Lambda = 0.01
    Model.n_ts_samples = 10
    Model.pts_of_int = numpy.concatenate([1, 4])
    Model.Guess_Q_C, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: q_matrix_comparisons(x, Model, Global), Model.Guess_Q_C, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    ##
    Model.ShowComparison = 'FullSnapshot'
    # Model.Guess_Q_C(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
    Outputs.Q.mean_err_C, Grad_calc, Outputs.Q.err_mat_C = q_matrix_comparisons(Model.Guess_Q_C, Model, Global,
                                                                                nargout=3)
    # print('=========================================================')
    # print('Composite:')
    # print(Outputs.Q.err_mat_C)
    # print(strcat('mean error:', num2str(mean(Outputs.Q.err_mat_C))))

    ##
    Model.Guess_P_FS = qvec_to_pmat(Model.Guess_Q_FS, Model)
    Model.Guess_P_FL = qvec_to_pmat(Model.Guess_Q_FL, Model)
    Model.Guess_P_TS = qvec_to_pmat(Model.Guess_Q_TS, Model)
    Model.Guess_P_C = qvec_to_pmat(Model.Guess_Q_C, Model)
    ##
    P_Lb = numpy.zeros(Model.num_basis)
    P_Lb[1 == numpy.eye(P_Lb.size())] = -math.inf
    P_Lb[:, -1 + 1] = 0
    P_beq = sparse(Model.num_basis + 1, 1)
    P_beq[-1] = 1
    P_Aeq = sparse(Model.num_basis + 1, Model.num_basis ** 2)
    for kk in range(1, (Model.num_basis + 1)).reshape(-1):
        P_Aeq[kk, range(((kk - 1) * Model.num_basis) + 1, (kk * Model.num_basis))] = 1

    ##

    test_init_Q_vec = numpy.zeros(Model.num_basis, Model.num_basis + 1)
    test_init_Q_vec[:, -1] = Model.Coeff_vecs[1]
    ##
    Model.ShowComparison = ('FullSnapshot')
    Model.Guess_P_FS, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), test_init_Q_vec, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    Model.Guess_P_FS, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), Model.Guess_P_FS, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    Model.Guess_P_FS, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), Model.Guess_P_FS, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    Model.ShowComparison = ('FullSnapshot')
    # Model.Guess_P_FS(:,end) = Model.Coeff_vecs{1};
    Outputs.P.mean_err_FS, Grad_calc, Outputs.P.err_mat_FS = p_matrix_comparisons(Model.Guess_P_FS, Model,
                                                                                  Global, nargout=3)
    # print('=========================================================')
    # print('Best we can hope for:')
    # print(Outputs.P.err_mat_FS)
    # print(strcat('mean error:', num2str(mean(Outputs.P.err_mat_FS))))
    ##
    Model.ShowComparison = ('FirstLastPoint')
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    Model.Guess_P_FL, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), test_init_Q_vec, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    Model.Guess_P_FL, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), Model.Guess_P_FL, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    Model.Guess_P_FL, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), Model.Guess_P_FL, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    Model.ShowComparison = 'FullSnapshot'
    # Model.Guess_P_FL(:,end) = Model.Coeff_vecs{1};
    Outputs.P.mean_err_FL, Grad_calc, Outputs.P.err_mat_FL = p_matrix_comparisons(Model.Guess_P_FL, Model,
                                                                                  Global, nargout=3)
    # print('=========================================================')
    # print('Worst case scenario:')
    # print(Outputs.P.err_mat_FL)
    # print(strcat('mean error:', num2str(mean(Outputs.P.err_mat_FL))))
    ##
    Model.ShowComparison = 'FullTimeSeries'
    Model.n_ts_samples = 100
    Model.Guess_P_TS, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), test_init_Q_vec, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    Model.Guess_P_TS, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), Model.Guess_P_TS, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    Model.Guess_P_TS, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), Model.Guess_P_TS, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    Model.ShowComparison = 'FullSnapshot'
    Model.Guess_P_TS[:, -1] = Model.Coeff_vecs[1]
    Outputs.P.mean_err_TS, Grad_calc, Outputs.P.err_mat_TS = p_matrix_comparisons(Model.Guess_P_TS, Model,
                                                                                  Global, nargout=3)
    # print('=========================================================')
    # print('Full Time Series:')
    # print(Outputs.P.err_mat_TS)
    # print(strcat('mean error:', num2str(mean(Outputs.P.err_mat_TS))))
    ##
    Model.ShowComparison = 'Composite'
    Model.Lambda = 0
    Model.n_ts_samples = 10
    Model.pts_of_int = numpy.concatenate([1, 4])
    Model.Guess_P_C, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), test_init_Q_vec, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    ##
    Model.ShowComparison = 'Composite'
    Model.Lambda = 0.01
    Model.n_ts_samples = 10
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    Model.Guess_P_C, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), Model.Guess_P_C, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    ##
    Model.ShowComparison = 'Composite'
    Model.Lambda = 0.01
    Model.n_ts_samples = 10
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    Model.Guess_P_C, _, tmp_exitflag, _, tmp_output = fmincon(
        lambda x=None: p_matrix_comparisons(x, Model, Global), Model.Guess_P_C, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    ##
    Model.ShowComparison = 'FullSnapshot'
    Outputs.P.mean_err_C, Grad_calc, Outputs.P.err_mat_C = p_matrix_comparisons(Model.Guess_P_C, Model, Global,
                                                                                nargout=3)
    # Checking Derivatives Q then P

    test_init_P_mat = qvec_to_pmat(test_init_Q_vec, Model)
    ##

    Model.ShowComparison = 'FullSnapshot'
    test_init_Q_vec = process_q_matrix(Model, test_init_Q_vec, Global)

    Model.ShowComparison = 'FullSnapshot'
    test_init_P_mat = process_p_matrix(Model, test_init_P_mat, Global)

    Model.ShowComparison = 'FirstLastPoint'
    test_init_Q_vec = process_q_matrix(Model, test_init_Q_vec, Global)

    Model.ShowComparison = 'FirstLastPoint'
    test_init_P_mat = process_q_matrix(Model, test_init_P_mat, Global)

    Model.ShowComparison = 'FullTimeSeries'
    Model.n_ts_samples = 100
    test_init_Q_vec = process_q_matrix(Model, test_init_Q_vec, Global)

    Model.ShowComparison = 'FullTimeSeries'
    Model.n_ts_samples = 100
    test_init_P_mat = process_p_matrix(Model, test_init_P_mat, Global)

    Model.ShowComparison = 'Composite'
    Model.Lambda = 0.25
    Model.n_ts_samples = 10
    test_init_Q_vec = process_q_matrix(Model, test_init_Q_vec, Global)

    Model.ShowComparison = 'Composite'
    Model.Lambda = 0.25
    Model.n_ts_samples = 10
    test_init_P_mat = process_p_matrix(Model, test_init_P_mat, Global)
