import numpy
import math
from scipy.io import loadmat
from sklearn.cluster._kmeans import k_means

if __name__ == "__main_":
    # Load processed data
    Model = {}
    Global = loadmat('SDE_ProcessedData.mat')
    # Rescale time to (0,1) interval
    Global.time_pts = (Global.time_pts / (Global.time_pts(-1) - Global.time_pts(1)))
    Model.norm_fac = 'Relative'
    # Choose number and locations of basis function
    Model.num_basis = 30
    __, Model.Basis_Loc = k_means(Global.all_pts_MixModel, Model.num_basis, 'Replicates', 5, nargout=2)
    Model.subfunc = {}
    # Define form of basis functions used.
    Model.subfunc.An = lambda n=None, eps=None: 2 * math.pi ** (n / 2) / gamma(n / 2)
    Model.subfunc.Vn = lambda n=None, eps=None: (eps ** n) * (math.pi ** (n / 2)) / gamma((n / 2) + 1)
    Model.subfunc.Intn = lambda n=None, eps=None: eps ** n / (n * (n + 1))
    Model.basis_nd = lambda x=None, n=None, eps=None: max(0, 1 - (vecnorm(x.T, 2) / eps)) / (
        (Model.subfunc.An(n, eps) * Model.subfunc.Intn(n, eps)))
    # Choose scale parameter for each basis function. Connect to (minimum) N nearest basis functions.
    Model.subfunc.connections = (2)

    Model.Basis_eps = (numpy.zeros(Model.num_basis, 1))
    for ii in range(1, Model.num_basis.reshape(-1).size()):
        print(ii)
        d_tmp = pdist2(Model.Basis_Loc[ii, :], Model.Basis_Loc)
        d_tmp = sort(d_tmp[d_tmp > 0])
        Model.Basis_eps[ii] = 1.00001 * d_tmp[Model.subfunc.connections]

    # Use Halton Node Placing for Monte Carlo integration
    p = haltonset(Global.Data_Dim)
    X0 = net(p, 10 ** 5)
    for ii in range(1, Model.num_basis.reshape(-1).size()):
        # print(ii)
        for jj in range(ii, Model.num_basis).reshape(-1):
            tmp = Package_MC_Integration(ii, jj, Model, Global, X0)
            Model.MassMat[ii, jj] = tmp
            Model.MassMat[jj, ii] = tmp

    Model.MassMat = (sparse(Model.MassMat))
    # Plot Mass Matrix and Inverse
    Model.InvMassMat = (inv(Model.MassMat))
    Model.InvMassMat = ((Model.InvMassMat + Model.InvMassMat.T) / 2)
    # Find coefficient vectors for each time pt.
    Model.Coeff_vecs = (cell(numel(Global.time_pts), 1))
    for vv in range(1, numel(Global.time_pts)).reshape(-1):
        Model.Coeff_vecs[vv] = numpy.zeros(Model.num_basis, 1)
        for ii in range(1, Model.num_basis).reshape(-1):
            tmp = Model.basis_nd(Global.sample_pts[vv] - Model.Basis_Loc(ii, :), Global.Data_Dim,
                                 Model.Basis_eps(ii))
            Model.Coeff_vecs[vv][ii] = sum(tmp)
        Model.Coeff_vecs[vv] = Model.Coeff_vecs[vv] / sum(Model.Coeff_vecs[vv])

    #
    Model.TimeSeries = (cell(numel(Global.time_pts), Global.Num_pts_each_time(1)))
    for ww in range(1, Global.Num_pts_each_time(1)).reshape(-1):
        print(ww)
        for vv in range(1, numel(Global.time_pts)).reshape(-1):
            Model.TimeSeries[vv, ww] = numpy.zeros(Model.num_basis, 1)
            for ii in range(1, Model.num_basis).reshape(-1):
                tmp = Model.basis_nd(Global.sample_pts[vv](ww, :) - Model.Basis_Loc(ii, :),
                                     Global.Data_Dim, Model.Basis_eps(ii))
                Model.TimeSeries[vv, ww][ii] = tmp
            Model.TimeSeries[vv, ww] = Model.TimeSeries[vv, ww] / sum(Model.TimeSeries[vv, ww])

    # Plot graphs of a few time points.
    G_tmp = graph(Model.MassMat)
    # ddd/Package_RunMe_Comparison.m:82
    G_tmp = rmedge(G_tmp, range(1, numnodes(G_tmp)), range(1, numnodes(G_tmp)))
    # ddd/Package_RunMe_Comparison.m:83
    for vv in numpy.concatenate([1, 9, 11]).reshape(-1):
        figure(vv)
        plot(G_tmp, 'XData', Model.Basis_Loc(:, 1), 'YData', Model.Basis_Loc(:, 2))
        hold('on')
        scatter(Model.Basis_Loc(:, 1), Model.Basis_Loc(:, 2), 50, Model.Coeff_vecs[vv], 'filled')
        axis(numpy.concatenate([0, 5, 0, 5]))

    # Plot graphs of a few time points.
    G_tmp = graph(Model.MassMat)
    G_tmp = rmedge(G_tmp, range(1, numnodes(G_tmp)), range(1, numnodes(G_tmp)))
    for vv in range(1, 11).reshape(-1):
        figure(vv)
        plot(G_tmp, 'XData', Model.Basis_Loc(:, 1), 'YData', Model.Basis_Loc(:, 2))
        hold('on')
        scatter(Model.Basis_Loc(:, 1), Model.Basis_Loc(:, 2), 50, Model.TimeSeries[vv, 4], 'filled')
        axis(numpy.concatenate([0, 5, 0, 5]))

    # Plot coefficient matrices
    CoeffMat = cell2mat(Model.Coeff_vecs.T)
    figure()
    imagesc(CoeffMat)
    drawnow
    #  Define useful quantities for later.
    Model.InvSums = (sum(Model.InvMassMat))
    Model.Indc = (sparse(Model.MassMat != 0))
    Indc_sum = sum(Model.Indc)
    Indc_zero = sparse(Model.MassMat == 0)
    Indc_vec = find(Model.MassMat != 0)
    Indc_diag = find(ismember(Indc_vec.T, range(1, Model.num_basis ** 2, (Model.num_basis + 1))))
    # Define initial guess at P matrix
    init_P = double(Model.Indc)
    for vv in range(1, Model.num_basis).reshape(-1):
        init_P[vv, vv] = 0
        init_P[vv, vv] = - sum(init_P(:, vv))

    clear('vv')
    # Calculate corresponding Q and vectorise non-zero entries
    init_Q = dot(Model.MassMat, init_P)
    init_Q[Indc_zero] = 0
    init_Q_vec = init_Q(logical_not(Indc_zero))
    # Specify effective dimension of problem.
    eff_var = numel(init_Q_vec)
    print(strcat('Effective number of parameters:', num2str(eff_var)))
    print(strcat('Eff. basis functions for full DDD:', num2str(sqrt(eff_var))))
    # Append initial condition onto end of vector.
    init_Q_vec[range(-1 + 1, -1 + Model.num_basis)] = Model.Coeff_vecs[1]
    # Creating constraints matrix for Aeq * init_Q_vec = beq
    Q_Aeq = sparse(Model.num_basis + 1, Model.num_basis + eff_var)
    run_tot_tmp = 1
    for kk in range(1, Model.num_basis).reshape(-1):
        ind_tmp = range(run_tot_tmp, run_tot_tmp + Indc_sum(kk) - 1)
        Q_Aeq[kk, ind_tmp] = Model.InvSums(Model.Indc(kk, :))
        run_tot_tmp = run_tot_tmp + Indc_sum(kk)

    clear('kk', 'ind_tmp', 'run_tot_tmp')
    Q_Aeq[Model.num_basis + 1, range((-1 - Model.num_basis + 1), -1)] = 1
    Q_beq = sparse(Model.num_basis + 1, 1)
    Q_beq[-1] = 1
    imagesc(numpy.concatenate(
        [numpy.concatenate([[sum(dot(Model.InvMassMat, init_Q)).T], [nan]]), dot(Q_Aeq, init_Q_vec), Q_beq]))
    # Create positivity constraint Alt * init_Q_vec <= blt
    Q_Alt = sparse(Model.num_basis + eff_var, Model.num_basis + eff_var)
    run_tot_tmp = 1
    for kk in range(1, Model.num_basis).reshape(-1):
        xx = find(Model.Indc(kk, :))
        yy = Model.InvMassMat(xx, xx)
        ind_tmp = range(run_tot_tmp, run_tot_tmp + Indc_sum(kk) - 1)
        Q_Alt[ind_tmp, ind_tmp] = - yy
        run_tot_tmp = run_tot_tmp + Indc_sum(kk)

    ind_tmp = range(run_tot_tmp, run_tot_tmp + Model.num_basis - 1)
    Q_Alt[ind_tmp, ind_tmp] = - eye(Model.num_basis)
    for kk in Indc_diag.reshape(-1):
        Q_Alt[kk, :] = - Q_Alt(kk, :)

    Q_blt = sparse(Model.num_basis + eff_var, 1)
    clear('kk', 'ind_tmp')
    ##

    opts = optimoptions('fmincon', 'Display', 'iter-detailed')
    opts.UseParallel = (true(1, 1))
    opts.MaxIterations = (3000)
    opts.SpecifyObjectiveGradient = (true(1, 1))
    opts.Algorithm = ('sqp')
    opts.MaxIterations = (300)
    test_init_Q_vec = (init_Q_vec)
    test_init_Q_vec[range(1, eff_var)] = 0
    Model.ShowComparison = ('FullSnapshot')
    Model.Guess_Q_FS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), test_init_Q_vec, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.Guess_Q_FS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), Model.Guess_Q_FS, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.Guess_Q_FS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), Model.Guess_Q_FS, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    Model.ShowComparison = ('FullSnapshot')
    # Model.Guess_Q_FS(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
    Outputs.Q.mean_err_FS, Grad_calc, Outputs.Q.err_mat_FS = Package_Q_Matrix_Comparisons(Model.Guess_Q_FS, Model,
                                                                                          Global, nargout=3)
    print('=========================================================')
    print('Best we can hope for:')
    print(Outputs.Q.err_mat_FS)
    print(strcat('mean error:', num2str(mean(Outputs.Q.err_mat_FS))))
    ##
    Model.ShowComparison = ('FirstLastPoint')
    # ddd/Package_RunMe_Comparison.m:194
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    # ddd/Package_RunMe_Comparison.m:195
    Model.Guess_Q_FL, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), test_init_Q_vec, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:196
    Model.Guess_Q_FL, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), Model.Guess_Q_FL, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:197
    Model.Guess_Q_FL, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), Model.Guess_Q_FL, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:198
    Model.ShowComparison = ('FullSnapshot')
    # ddd/Package_RunMe_Comparison.m:199
    # Model.Guess_Q_FL(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
    Outputs.Q.mean_err_FL, Grad_calc, Outputs.Q.err_mat_FL = Package_Q_Matrix_Comparisons(Model.Guess_Q_FL, Model,
                                                                                          Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:201
    print('=========================================================')
    print('Worst case scenario:')
    print(Outputs.Q.err_mat_FL)
    print(strcat('mean error:', num2str(mean(Outputs.Q.err_mat_FL))))
    ##
    Model.ShowComparison = ('FullTimeSeries')
    # ddd/Package_RunMe_Comparison.m:207
    Model.n_ts_samples = (100)
    # ddd/Package_RunMe_Comparison.m:208
    Model.Guess_Q_TS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), test_init_Q_vec, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:209
    Model.Guess_Q_TS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), Model.Guess_Q_TS, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:210
    Model.Guess_Q_TS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), Model.Guess_Q_TS, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:211
    Model.ShowComparison = ('FullSnapshot')
    # ddd/Package_RunMe_Comparison.m:212
    Model.Guess_Q_TS[range(-1 - Model.num_basis + 1, -1)] = Model.Coeff_vecs[1]
    # ddd/Package_RunMe_Comparison.m:213
    Outputs.Q.mean_err_TS, Grad_calc, Outputs.Q.err_mat_TS = Package_Q_Matrix_Comparisons(Model.Guess_Q_TS, Model,
                                                                                          Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:214
    print('=========================================================')
    print('Full Time Series:')
    print(Outputs.Q.err_mat_TS)
    print(strcat('mean error:', num2str(mean(Outputs.Q.err_mat_TS))))
    ##
    Model.ShowComparison = ('Composite')
    # ddd/Package_RunMe_Comparison.m:221
    Model.Lambda = (0)
    # ddd/Package_RunMe_Comparison.m:222
    Model.n_ts_samples = (10)
    # ddd/Package_RunMe_Comparison.m:222
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    # ddd/Package_RunMe_Comparison.m:222
    Model.Guess_Q_C, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), test_init_Q_vec, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:223
    ##
    Model.ShowComparison = ('Composite')
    # ddd/Package_RunMe_Comparison.m:225
    Model.Lambda = (0.01)
    # ddd/Package_RunMe_Comparison.m:226
    Model.n_ts_samples = (10)
    # ddd/Package_RunMe_Comparison.m:226
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    # ddd/Package_RunMe_Comparison.m:226
    Model.Guess_Q_C, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), Model.Guess_Q_C, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:227
    ##
    Model.ShowComparison = ('Composite')
    # ddd/Package_RunMe_Comparison.m:229
    Model.Lambda = (0.01)
    # ddd/Package_RunMe_Comparison.m:230
    Model.n_ts_samples = (10)
    # ddd/Package_RunMe_Comparison.m:230
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    # ddd/Package_RunMe_Comparison.m:230
    Model.Guess_Q_C, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_Q_Matrix_Comparisons(x, Model, Global), Model.Guess_Q_C, Q_Alt, Q_blt, Q_Aeq, Q_beq, [],
        [], [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:231
    ##
    Model.ShowComparison = ('FullSnapshot')
    # ddd/Package_RunMe_Comparison.m:233
    # Model.Guess_Q_C(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
    Outputs.Q.mean_err_C, Grad_calc, Outputs.Q.err_mat_C = Package_Q_Matrix_Comparisons(Model.Guess_Q_C, Model, Global,
                                                                                        nargout=3)
    # ddd/Package_RunMe_Comparison.m:235
    print('=========================================================')
    print('Composite:')
    print(Outputs.Q.err_mat_C)
    print(strcat('mean error:', num2str(mean(Outputs.Q.err_mat_C))))

    ##
    figure(1)
    plot(log10(1 + dot(89, Global.time_pts)), 2 + log10(sqrt(Outputs.Q.err_mat_TS)), '--x', 'LineWidth', 2)
    hold('on')
    plot(log10(1 + dot(89, Global.time_pts)), 2 + log10(sqrt(Outputs.Q.err_mat_FS)), '--x', 'LineWidth', 2)
    hold('on')
    plot(log10(1 + dot(89, Global.time_pts)), 2 + log10(sqrt(Outputs.Q.err_mat_FL)), '--x', 'LineWidth', 2)
    hold('on')
    plot(log10(1 + dot(89, Global.time_pts)), 2 + log10(sqrt(Outputs.Q.err_mat_C)), '--x', 'LineWidth', 2)
    hold('on')
    axis(numpy.concatenate([- 0.25, 2, - 0.5, 2.5]))
    legend('Trajectory', 'Full Snapshot', 'Two Snapshot', 'Integrated')
    # set(gca,'XTick',0:2)
    # set(gca,'YTick',-3:3)
    set(gca, 'fontsize', 17)
    axis('square')
    tightfig()
    ##
    Model.Guess_P_FS = (Package_Func_QVec2PMat(Model.Guess_Q_FS, Model))
    # ddd/Package_RunMe_Comparison.m:255
    Model.Guess_P_FL = (Package_Func_QVec2PMat(Model.Guess_Q_FL, Model))
    # ddd/Package_RunMe_Comparison.m:256
    Model.Guess_P_TS = (Package_Func_QVec2PMat(Model.Guess_Q_TS, Model))
    # ddd/Package_RunMe_Comparison.m:257
    Model.Guess_P_C = (Package_Func_QVec2PMat(Model.Guess_Q_C, Model))
    # ddd/Package_RunMe_Comparison.m:258
    ##
    P_Lb = numpy.zeros(Model.num_basis)
    # ddd/Package_RunMe_Comparison.m:261
    P_Lb[1 == eye(size(P_Lb))] = - Inf
    # ddd/Package_RunMe_Comparison.m:262
    P_Lb[:, -1 + 1] = 0
    # ddd/Package_RunMe_Comparison.m:263
    P_beq = sparse(Model.num_basis + 1, 1)
    # ddd/Package_RunMe_Comparison.m:264
    P_beq[-1] = 1
    # ddd/Package_RunMe_Comparison.m:265
    P_Aeq = sparse(Model.num_basis + 1, Model.num_basis ** 2)
    # ddd/Package_RunMe_Comparison.m:266
    for kk in range(1, (Model.num_basis + 1)).reshape(-1):
        P_Aeq[kk, range(((kk - 1) * Model.num_basis) + 1, (kk * Model.num_basis))] = 1
    # ddd/Package_RunMe_Comparison.m:268

    ##

    test_init_Q_vec = numpy.zeros(Model.num_basis, Model.num_basis + 1)
    # ddd/Package_RunMe_Comparison.m:272
    test_init_Q_vec[:, -1] = Model.Coeff_vecs[1]
    # ddd/Package_RunMe_Comparison.m:273
    ##
    Model.ShowComparison = ('FullSnapshot')
    # ddd/Package_RunMe_Comparison.m:276
    Model.Guess_P_FS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), test_init_Q_vec, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:277
    Model.Guess_P_FS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), Model.Guess_P_FS, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:278
    Model.Guess_P_FS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), Model.Guess_P_FS, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:279
    Model.ShowComparison = ('FullSnapshot')
    # ddd/Package_RunMe_Comparison.m:280
    # Model.Guess_P_FS(:,end) = Model.Coeff_vecs{1};
    Outputs.P.mean_err_FS, Grad_calc, Outputs.P.err_mat_FS = Package_P_Matrix_Comparisons(Model.Guess_P_FS, Model,
                                                                                          Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:282
    print('=========================================================')
    print('Best we can hope for:')
    print(Outputs.P.err_mat_FS)
    print(strcat('mean error:', num2str(mean(Outputs.P.err_mat_FS))))
    ##
    Model.ShowComparison = ('FirstLastPoint')
    # ddd/Package_RunMe_Comparison.m:288
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    # ddd/Package_RunMe_Comparison.m:289
    Model.Guess_P_FL, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), test_init_Q_vec, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:290
    Model.Guess_P_FL, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), Model.Guess_P_FL, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:291
    Model.Guess_P_FL, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), Model.Guess_P_FL, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:292
    Model.ShowComparison = ('FullSnapshot')
    # ddd/Package_RunMe_Comparison.m:293
    # Model.Guess_P_FL(:,end) = Model.Coeff_vecs{1};
    Outputs.P.mean_err_FL, Grad_calc, Outputs.P.err_mat_FL = Package_P_Matrix_Comparisons(Model.Guess_P_FL, Model,
                                                                                          Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:295
    print('=========================================================')
    print('Worst case scenario:')
    print(Outputs.P.err_mat_FL)
    print(strcat('mean error:', num2str(mean(Outputs.P.err_mat_FL))))
    ##
    Model.ShowComparison = ('FullTimeSeries')
    # ddd/Package_RunMe_Comparison.m:301
    Model.n_ts_samples = (100)
    # ddd/Package_RunMe_Comparison.m:302
    Model.Guess_P_TS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), test_init_Q_vec, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:303
    Model.Guess_P_TS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), Model.Guess_P_TS, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:304
    Model.Guess_P_TS, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), Model.Guess_P_TS, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:305
    Model.ShowComparison = ('FullSnapshot')
    # ddd/Package_RunMe_Comparison.m:306
    Model.Guess_P_TS[:, -1] = Model.Coeff_vecs[1]
    # ddd/Package_RunMe_Comparison.m:307
    Outputs.P.mean_err_TS, Grad_calc, Outputs.P.err_mat_TS = Package_P_Matrix_Comparisons(Model.Guess_P_TS, Model,
                                                                                          Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:308
    print('=========================================================')
    print('Full Time Series:')
    print(Outputs.P.err_mat_TS)
    print(strcat('mean error:', num2str(mean(Outputs.P.err_mat_TS))))
    ##
    Model.ShowComparison = ('Composite')
    # ddd/Package_RunMe_Comparison.m:315
    Model.Lambda = (0)
    # ddd/Package_RunMe_Comparison.m:316
    Model.n_ts_samples = (10)
    # ddd/Package_RunMe_Comparison.m:316
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    # ddd/Package_RunMe_Comparison.m:316
    Model.Guess_P_C, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), test_init_Q_vec, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:317
    ##
    Model.ShowComparison = ('Composite')
    # ddd/Package_RunMe_Comparison.m:319
    Model.Lambda = (0.01)
    # ddd/Package_RunMe_Comparison.m:320
    Model.n_ts_samples = (10)
    # ddd/Package_RunMe_Comparison.m:320
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    # ddd/Package_RunMe_Comparison.m:320
    Model.Guess_P_C, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), Model.Guess_P_C, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:321
    ##
    Model.ShowComparison = ('Composite')
    # ddd/Package_RunMe_Comparison.m:323
    Model.Lambda = (0.01)
    # ddd/Package_RunMe_Comparison.m:324
    Model.n_ts_samples = (10)
    # ddd/Package_RunMe_Comparison.m:324
    Model.pts_of_int = (numpy.concatenate([1, 4]))
    # ddd/Package_RunMe_Comparison.m:324
    Model.Guess_P_C, __, tmp_exitflag, __, tmp_output = fmincon(
        lambda x=None: Package_P_Matrix_Comparisons(x, Model, Global), Model.Guess_P_C, [], [], P_Aeq, P_beq, P_Lb, [],
        [], opts, nargout=5)
    # ddd/Package_RunMe_Comparison.m:325
    ##
    Model.ShowComparison = ('FullSnapshot')
    # ddd/Package_RunMe_Comparison.m:327
    # Model.Guess_P_C(:,end) = Model.Coeff_vecs{1};
    Outputs.P.mean_err_C, Grad_calc, Outputs.P.err_mat_C = Package_P_Matrix_Comparisons(Model.Guess_P_C, Model, Global,
                                                                                        nargout=3)
    # ddd/Package_RunMe_Comparison.m:329
    print('=========================================================')
    print('Composite:')
    print(Outputs.P.err_mat_C)
    print(strcat('mean error:', num2str(mean(Outputs.P.err_mat_C))))
    ##
    figure(2)
    plot(log10(1 + dot(89, Global.time_pts)), 2 + log10(sqrt(Outputs.P.err_mat_TS)), '--x', 'LineWidth', 2)
    hold('on')
    plot(log10(1 + dot(89, Global.time_pts)), 2 + log10(sqrt(Outputs.P.err_mat_FS)), '--x', 'LineWidth', 2)
    hold('on')
    plot(log10(1 + dot(89, Global.time_pts)), 2 + log10(sqrt(Outputs.P.err_mat_FL)), '--x', 'LineWidth', 2)
    hold('on')
    plot(log10(1 + dot(89, Global.time_pts)), 2 + log10(sqrt(Outputs.P.err_mat_C)), '--x', 'LineWidth', 2)
    hold('on')
    axis(numpy.concatenate([- 0.25, 2, - 0.5, 2.5]))
    legend('Trajectory', 'Full Snapshot', 'Two Snapshot', 'Integrated')
    # set(gca,'XTick',0:2)
    # set(gca,'YTick',-3:3)
    set(gca, 'fontsize', 17)
    axis('square')
    tightfig()
    ##

    ##
    print('All Time Series')
    print(mean(sqrt(Outputs.Q.err_mat_TS(range(2, -1)))))
    print(mean(sqrt(Outputs.P.err_mat_TS(range(2, -1)))))
    ##
    print('All Snapshots')
    print(mean(sqrt(Outputs.Q.err_mat_FS)))
    print(mean(sqrt(Outputs.P.err_mat_FS)))
    ##
    print('Two Snapshots')
    print(mean(sqrt(Outputs.Q.err_mat_FL)))
    print(mean(sqrt(Outputs.P.err_mat_FL)))
    ##
    print('Composite')
    print(mean(sqrt(Outputs.Q.err_mat_C)))
    print(mean(sqrt(Outputs.P.err_mat_C)))
    ##

    ##
    ####### Checking Derivatives Q then P

    test_init_P_mat = Package_Func_QVec2PMat(test_init_Q_vec, Model)
    # ddd/Package_RunMe_Comparison.m:378
    ##

    format('long')
    Model.ShowComparison = ('FullSnapshot')
    # ddd/Package_RunMe_Comparison.m:382
    a, b, c = Package_Q_Matrix_Comparisons(test_init_Q_vec, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:383
    test_init_Q_vec[2] = test_init_Q_vec(2) + 0.0001
    # ddd/Package_RunMe_Comparison.m:384
    A, B, C = Package_Q_Matrix_Comparisons(test_init_Q_vec, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:385
    print(strcat('Derivative checking:', num2str(numpy.concatenate([(A - a) / 0.0001, b(2), B(2)]))))
    print(strcat('Error:', num2str(a)))
    format('short')

    format('long')
    Model.ShowComparison = ('FullSnapshot')
    # ddd/Package_RunMe_Comparison.m:391
    a, b, __ = Package_P_Matrix_Comparisons(test_init_P_mat, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:392
    test_init_P_mat[2] = test_init_P_mat(2) + 0.0001
    # ddd/Package_RunMe_Comparison.m:393
    A, B, __ = Package_P_Matrix_Comparisons(test_init_P_mat, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:394
    print(strcat('Derivative checking:', num2str(numpy.concatenate([(A - a) / 0.0001, b(2), B(2)]))))
    print(strcat('Error:', num2str(a)))
    format('short')
    ##

    format('long')
    Model.ShowComparison = ('FirstLastPoint')
    # ddd/Package_RunMe_Comparison.m:401
    Model.pts_of_int = (numpy.concatenate([1, 6]))
    # ddd/Package_RunMe_Comparison.m:402
    a, b, __ = Package_Q_Matrix_Comparisons(test_init_Q_vec, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:403
    test_init_Q_vec[2] = test_init_Q_vec(2) + 0.0001
    # ddd/Package_RunMe_Comparison.m:404
    A, B, __ = Package_Q_Matrix_Comparisons(test_init_Q_vec, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:405
    print(strcat('Derivative checking:', num2str(numpy.concatenate([(A - a) / 0.0001, b(2), B(2)]))))
    print(strcat('Error:', num2str(a)))
    format('short')

    format('long')
    Model.ShowComparison = ('FirstLastPoint')
    # ddd/Package_RunMe_Comparison.m:411
    Model.pts_of_int = (numpy.concatenate([1, 6]))
    # ddd/Package_RunMe_Comparison.m:412
    a, b, __ = Package_P_Matrix_Comparisons(test_init_P_mat, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:413
    test_init_P_mat[2] = test_init_P_mat(2) + 0.0001
    # ddd/Package_RunMe_Comparison.m:414
    A, B, __ = Package_P_Matrix_Comparisons(test_init_P_mat, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:415
    print(strcat('Derivative checking:', num2str(numpy.concatenate([(A - a) / 0.0001, b(2), B(2)]))))
    print(strcat('Error:', num2str(a)))
    format('short')
    ##

    format('long')
    Model.ShowComparison = ('FullTimeSeries')
    # ddd/Package_RunMe_Comparison.m:422
    Model.n_ts_samples = (100)
    # ddd/Package_RunMe_Comparison.m:423
    a, b, __ = Package_Q_Matrix_Comparisons(test_init_Q_vec, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:424
    test_init_Q_vec[2] = test_init_Q_vec(2) + 0.0001
    # ddd/Package_RunMe_Comparison.m:425
    A, B, __ = Package_Q_Matrix_Comparisons(test_init_Q_vec, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:426
    print(strcat('Derivative checking:', num2str(numpy.concatenate([(A - a) / 0.0001, b(2), B(2)]))))
    print(strcat('Error:', num2str(a)))
    format('short')

    format('long')
    Model.ShowComparison = ('FullTimeSeries')
    # ddd/Package_RunMe_Comparison.m:432
    Model.n_ts_samples = (100)
    # ddd/Package_RunMe_Comparison.m:433
    a, b, __ = Package_P_Matrix_Comparisons(test_init_P_mat, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:434
    test_init_P_mat[2] = test_init_P_mat(2) + 0.0001
    # ddd/Package_RunMe_Comparison.m:435
    A, B, __ = Package_P_Matrix_Comparisons(test_init_P_mat, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:436
    print(strcat('Derivative checking:', num2str(numpy.concatenate([(A - a) / 0.0001, b(2), B(2)]))))
    print(strcat('Error:', num2str(a)))
    format('short')
    ##

    format('long')
    Model.ShowComparison = ('Composite')
    # ddd/Package_RunMe_Comparison.m:443
    Model.Lambda = (0.25)
    # ddd/Package_RunMe_Comparison.m:444
    Model.n_ts_samples = (10)
    # ddd/Package_RunMe_Comparison.m:444
    Model.pts_of_int = (numpy.concatenate([1, 6]))
    # ddd/Package_RunMe_Comparison.m:444
    a, b, __ = Package_Q_Matrix_Comparisons(test_init_Q_vec, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:445
    test_init_Q_vec[2] = test_init_Q_vec(2) + 0.0001
    # ddd/Package_RunMe_Comparison.m:446
    A, B, __ = Package_Q_Matrix_Comparisons(test_init_Q_vec, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:447
    print(strcat('Derivative checking:', num2str(numpy.concatenate([(A - a) / 0.0001, b(2), B(2)]))))
    print(strcat('Error:', num2str(a)))
    format('short')

    format('long')
    Model.ShowComparison = ('Composite')
    # ddd/Package_RunMe_Comparison.m:453
    Model.Lambda = (0.25)
    # ddd/Package_RunMe_Comparison.m:454
    Model.n_ts_samples = (10)
    # ddd/Package_RunMe_Comparison.m:454
    Model.pts_of_int = (numpy.concatenate([1, 6]))
    # ddd/Package_RunMe_Comparison.m:454
    a, b, __ = Package_P_Matrix_Comparisons(test_init_P_mat, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:455
    test_init_P_mat[2] = test_init_P_mat(2) + 0.0001
    # ddd/Package_RunMe_Comparison.m:456
    A, B, __ = Package_P_Matrix_Comparisons(test_init_P_mat, Model, Global, nargout=3)
    # ddd/Package_RunMe_Comparison.m:457
    print(strcat('Derivative checking:', num2str(numpy.concatenate([(A - a) / 0.0001, b(2), B(2)]))))
    print(strcat('Error:', num2str(a)))
    format('short')
