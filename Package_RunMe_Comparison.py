# Generated with SMOP  0.41-beta
from libsmop import *
# ddd/Package_RunMe_Comparison.m
if __name__ =="__main_":
    # Load processed data
    Model = {}
    load('SDE_ProcessedData.mat')
    # Rescale time to (0,1) interval
    Global.time_pts = copy(Global.time_pts / (Global.time_pts(end()) - Global.time_pts(1)))
# ddd/Package_RunMe_Comparison.m:7
    Model.norm_fac = copy('Relative')
# ddd/Package_RunMe_Comparison.m:8
    ## Choose number and locations of basis function
    Model.num_basis = copy(30)
# ddd/Package_RunMe_Comparison.m:10
    __,Model.Basis_Loc=kmeans(Global.all_pts_MixModel,Model.num_basis,'Replicates',5,nargout=2)
# ddd/Package_RunMe_Comparison.m:11
    Model.subfunc ={}
    ## Define form of basis functions used.
    Model.subfunc.An = copy(lambda n=None,eps=None: dot(2,pi ** (n / 2)) / gamma(n / 2))
# ddd/Package_RunMe_Comparison.m:14
    Model.subfunc.Vn = copy(lambda n=None,eps=None: multiply((eps ** n),pi ** (n / 2)) / gamma((n / 2) + 1))
# ddd/Package_RunMe_Comparison.m:15
    Model.subfunc.Intn = copy(lambda n=None,eps=None: eps ** n / (dot(n,(n + 1))))
# ddd/Package_RunMe_Comparison.m:16
    Model.basis_nd = copy(lambda x=None,n=None,eps=None: max(0,1 - (vecnorm(x.T,2) / eps)) / (dot(Model.subfunc.An(n,eps),Model.subfunc.Intn(n,eps))))
# ddd/Package_RunMe_Comparison.m:17
    ## Choose scale parameter for each basis function. Connect to (minimum) N nearest basis functions.
    Model.subfunc.connections = copy(2)
# ddd/Package_RunMe_Comparison.m:20
    
    Model.Basis_eps = copy(zeros(Model.num_basis,1))
# ddd/Package_RunMe_Comparison.m:21
    for ii in arange(1,Model.num_basis).reshape(-1):
        disp(ii)
        d_tmp=pdist2(Model.Basis_Loc(ii,arange()),Model.Basis_Loc)
# ddd/Package_RunMe_Comparison.m:24
        d_tmp=sort(d_tmp(d_tmp > 0))
# ddd/Package_RunMe_Comparison.m:25
        Model.Basis_eps[ii]=dot(1.00001,d_tmp(Model.subfunc.connections))
# ddd/Package_RunMe_Comparison.m:26
    
    clear('d_tmp','ii')
    ## Use Halton Node Placing for Monte Carlo integration
    p=haltonset(Global.Data_Dim)
# ddd/Package_RunMe_Comparison.m:31
    X0=net(p,10 ** 5)
# ddd/Package_RunMe_Comparison.m:32
    for ii in arange(1,Model.num_basis).reshape(-1):
        disp(ii)
        for jj in arange(ii,Model.num_basis).reshape(-1):
            tmp=Package_MC_Integration(ii,jj,Model,Global,X0)
# ddd/Package_RunMe_Comparison.m:36
            Model.MassMat[ii,jj]=tmp
# ddd/Package_RunMe_Comparison.m:37
            Model.MassMat[jj,ii]=tmp
# ddd/Package_RunMe_Comparison.m:38
    
    clear('ii','jj','tmp','X0','p')
    Model.MassMat = copy(sparse(Model.MassMat))
# ddd/Package_RunMe_Comparison.m:42
    ## Plot Mass Matrix and Inverse
    figure()
    imagesc(log(Model.MassMat))
    Model.InvMassMat = copy(inv(Model.MassMat))
# ddd/Package_RunMe_Comparison.m:48
    Model.InvMassMat = copy((Model.InvMassMat + Model.InvMassMat.T) / 2)
# ddd/Package_RunMe_Comparison.m:49
    figure()
    imagesc(log(abs(Model.InvMassMat)))
    ## Find coefficient vectors for each time pt.
    Model.Coeff_vecs = copy(cell(numel(Global.time_pts),1))
# ddd/Package_RunMe_Comparison.m:55
    for vv in arange(1,numel(Global.time_pts)).reshape(-1):
        disp(vv)
        Model.Coeff_vecs[vv]=zeros(Model.num_basis,1)
# ddd/Package_RunMe_Comparison.m:58
        for ii in arange(1,Model.num_basis).reshape(-1):
            tmp=Model.basis_nd(Global.sample_pts[vv] - Model.Basis_Loc(ii,arange()),Global.Data_Dim,Model.Basis_eps(ii))
# ddd/Package_RunMe_Comparison.m:60
            Model.Coeff_vecs[vv][ii]=sum(tmp)
# ddd/Package_RunMe_Comparison.m:61
        Model.Coeff_vecs[vv]=Model.Coeff_vecs[vv] / sum(Model.Coeff_vecs[vv])
# ddd/Package_RunMe_Comparison.m:63
    
    clear('vv','ii','tmp')
    ##
    Model.TimeSeries = copy(cell(numel(Global.time_pts),Global.Num_pts_each_time(1)))
# ddd/Package_RunMe_Comparison.m:67
    for ww in arange(1,Global.Num_pts_each_time(1)).reshape(-1):
        disp(ww)
        for vv in arange(1,numel(Global.time_pts)).reshape(-1):
            #disp(vv)
            Model.TimeSeries[vv,ww]=zeros(Model.num_basis,1)
# ddd/Package_RunMe_Comparison.m:72
            for ii in arange(1,Model.num_basis).reshape(-1):
                tmp=Model.basis_nd(Global.sample_pts[vv](ww,arange()) - Model.Basis_Loc(ii,arange()),Global.Data_Dim,Model.Basis_eps(ii))
# ddd/Package_RunMe_Comparison.m:74
                Model.TimeSeries[vv,ww][ii]=tmp
# ddd/Package_RunMe_Comparison.m:75
            Model.TimeSeries[vv,ww]=Model.TimeSeries[vv,ww] / sum(Model.TimeSeries[vv,ww])
# ddd/Package_RunMe_Comparison.m:77
    
    clear('vv','ww','ii','tmp')
    ## Plot graphs of a few time points.
    G_tmp=graph(Model.MassMat)
# ddd/Package_RunMe_Comparison.m:82
    G_tmp=rmedge(G_tmp,arange(1,numnodes(G_tmp)),arange(1,numnodes(G_tmp)))
# ddd/Package_RunMe_Comparison.m:83
    for vv in concat([1,9,11]).reshape(-1):
        figure(vv)
        plot(G_tmp,'XData',Model.Basis_Loc(arange(),1),'YData',Model.Basis_Loc(arange(),2))
        hold('on')
        scatter(Model.Basis_Loc(arange(),1),Model.Basis_Loc(arange(),2),50,Model.Coeff_vecs[vv],'filled')
        axis(concat([0,5,0,5]))
    
    clear('G_tmp','vv')
    ## Plot graphs of a few time points.
    G_tmp=graph(Model.MassMat)
# ddd/Package_RunMe_Comparison.m:92
    G_tmp=rmedge(G_tmp,arange(1,numnodes(G_tmp)),arange(1,numnodes(G_tmp)))
# ddd/Package_RunMe_Comparison.m:93
    for vv in arange(1,11).reshape(-1):
        figure(vv)
        plot(G_tmp,'XData',Model.Basis_Loc(arange(),1),'YData',Model.Basis_Loc(arange(),2))
        hold('on')
        scatter(Model.Basis_Loc(arange(),1),Model.Basis_Loc(arange(),2),50,Model.TimeSeries[vv,4],'filled')
        axis(concat([0,5,0,5]))
    
    clear('G_tmp','vv')
    ## Plot coefficient matrices
    CoeffMat=cell2mat(Model.Coeff_vecs.T)
# ddd/Package_RunMe_Comparison.m:104
    figure()
    imagesc(CoeffMat)
    drawnow
    ##  Define useful quantities for later.
    Model.InvSums = copy(sum(Model.InvMassMat))
# ddd/Package_RunMe_Comparison.m:110
    Model.Indc = copy(sparse(Model.MassMat != 0))
# ddd/Package_RunMe_Comparison.m:111
    Indc_sum=sum(Model.Indc)
# ddd/Package_RunMe_Comparison.m:112
    Indc_zero=sparse(Model.MassMat == 0)
# ddd/Package_RunMe_Comparison.m:113
    Indc_vec=find(Model.MassMat != 0)
# ddd/Package_RunMe_Comparison.m:114
    Indc_diag=find(ismember(Indc_vec.T,arange(1,Model.num_basis ** 2,(Model.num_basis + 1))))
# ddd/Package_RunMe_Comparison.m:115
    ## Define initial guess at P matrix
    init_P=double(Model.Indc)
# ddd/Package_RunMe_Comparison.m:119
    for vv in arange(1,Model.num_basis).reshape(-1):
        init_P[vv,vv]=0
# ddd/Package_RunMe_Comparison.m:121
        init_P[vv,vv]=- sum(init_P(arange(),vv))
# ddd/Package_RunMe_Comparison.m:122
    
    clear('vv')
    ## Calculate corresponding Q and vectorise non-zero entries
    init_Q=dot(Model.MassMat,init_P)
# ddd/Package_RunMe_Comparison.m:126
    init_Q[Indc_zero]=0
# ddd/Package_RunMe_Comparison.m:127
    init_Q_vec=init_Q(logical_not(Indc_zero))
# ddd/Package_RunMe_Comparison.m:128
    ## Specify effective dimension of problem.
    eff_var=numel(init_Q_vec)
# ddd/Package_RunMe_Comparison.m:130
    disp(strcat('Effective number of parameters:',num2str(eff_var)))
    disp(strcat('Eff. basis functions for full DDD:',num2str(sqrt(eff_var))))
    ## Append initial condition onto end of vector.
    init_Q_vec[arange(end() + 1,end() + Model.num_basis)]=Model.Coeff_vecs[1]
# ddd/Package_RunMe_Comparison.m:134
    ## Creating constraints matrix for Aeq * init_Q_vec = beq
    Q_Aeq=sparse(Model.num_basis + 1,Model.num_basis + eff_var)
# ddd/Package_RunMe_Comparison.m:136
    run_tot_tmp=1
# ddd/Package_RunMe_Comparison.m:137
    for kk in arange(1,Model.num_basis).reshape(-1):
        ind_tmp=arange(run_tot_tmp,run_tot_tmp + Indc_sum(kk) - 1)
# ddd/Package_RunMe_Comparison.m:139
        Q_Aeq[kk,ind_tmp]=Model.InvSums(Model.Indc(kk,arange()))
# ddd/Package_RunMe_Comparison.m:140
        run_tot_tmp=run_tot_tmp + Indc_sum(kk)
# ddd/Package_RunMe_Comparison.m:141
    
    clear('kk','ind_tmp','run_tot_tmp')
    Q_Aeq[Model.num_basis + 1,arange((end() - Model.num_basis + 1),end())]=1
# ddd/Package_RunMe_Comparison.m:144
    Q_beq=sparse(Model.num_basis + 1,1)
# ddd/Package_RunMe_Comparison.m:145
    Q_beq[end()]=1
# ddd/Package_RunMe_Comparison.m:146
    ##
    imagesc(concat([concat([[sum(dot(Model.InvMassMat,init_Q)).T],[nan]]),dot(Q_Aeq,init_Q_vec),Q_beq]))
    ## Create positivity constraint Alt * init_Q_vec <= blt
    Q_Alt=sparse(Model.num_basis + eff_var,Model.num_basis + eff_var)
# ddd/Package_RunMe_Comparison.m:151
    run_tot_tmp=1
# ddd/Package_RunMe_Comparison.m:152
    for kk in arange(1,Model.num_basis).reshape(-1):
        xx=find(Model.Indc(kk,arange()))
# ddd/Package_RunMe_Comparison.m:154
        yy=Model.InvMassMat(xx,xx)
# ddd/Package_RunMe_Comparison.m:155
        ind_tmp=arange(run_tot_tmp,run_tot_tmp + Indc_sum(kk) - 1)
# ddd/Package_RunMe_Comparison.m:156
        Q_Alt[ind_tmp,ind_tmp]=- yy
# ddd/Package_RunMe_Comparison.m:157
        run_tot_tmp=run_tot_tmp + Indc_sum(kk)
# ddd/Package_RunMe_Comparison.m:158
    
    ind_tmp=arange(run_tot_tmp,run_tot_tmp + Model.num_basis - 1)
# ddd/Package_RunMe_Comparison.m:160
    Q_Alt[ind_tmp,ind_tmp]=- eye(Model.num_basis)
# ddd/Package_RunMe_Comparison.m:161
    for kk in Indc_diag.reshape(-1):
        Q_Alt[kk,arange()]=- Q_Alt(kk,arange())
# ddd/Package_RunMe_Comparison.m:163
    
    Q_blt=sparse(Model.num_basis + eff_var,1)
# ddd/Package_RunMe_Comparison.m:165
    clear('kk','ind_tmp')
    ##
    
    opts=optimoptions('fmincon','Display','iter-detailed')
# ddd/Package_RunMe_Comparison.m:169
    opts.UseParallel = copy(true(1,1))
# ddd/Package_RunMe_Comparison.m:170
    opts.MaxIterations = copy(3000)
# ddd/Package_RunMe_Comparison.m:171
    opts.SpecifyObjectiveGradient = copy(true(1,1))
# ddd/Package_RunMe_Comparison.m:172
    opts.Algorithm = copy('sqp')
# ddd/Package_RunMe_Comparison.m:173
    #opts.TypicalX = init_K;
    opts.MaxIterations = copy(300)
# ddd/Package_RunMe_Comparison.m:175
    ##
#[Model.Guess_Q,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Find_Q_Matrix(x,Model,Global),init_Q_vec,Alt,blt,Aeq,beq,[],[],[],opts);
##
    test_init_Q_vec=copy(init_Q_vec)
# ddd/Package_RunMe_Comparison.m:179
    test_init_Q_vec[arange(1,eff_var)]=0
# ddd/Package_RunMe_Comparison.m:180
    ##
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:182
    Model.Guess_Q_FS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:183
    Model.Guess_Q_FS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_FS,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:184
    Model.Guess_Q_FS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_FS,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:185
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:186
    #Model.Guess_Q_FS(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
    Outputs.Q.mean_err_FS,Grad_calc,Outputs.Q.err_mat_FS=Package_Q_Matrix_Comparisons(Model.Guess_Q_FS,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:188
    disp('=========================================================')
    disp('Best we can hope for:')
    disp(Outputs.Q.err_mat_FS)
    disp(strcat('mean error:',num2str(mean(Outputs.Q.err_mat_FS))))
    ##
    Model.ShowComparison = copy('FirstLastPoint')
# ddd/Package_RunMe_Comparison.m:194
    Model.pts_of_int = copy(concat([1,4]))
# ddd/Package_RunMe_Comparison.m:195
    Model.Guess_Q_FL,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:196
    Model.Guess_Q_FL,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_FL,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:197
    Model.Guess_Q_FL,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_FL,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:198
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:199
    #Model.Guess_Q_FL(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
    Outputs.Q.mean_err_FL,Grad_calc,Outputs.Q.err_mat_FL=Package_Q_Matrix_Comparisons(Model.Guess_Q_FL,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:201
    disp('=========================================================')
    disp('Worst case scenario:')
    disp(Outputs.Q.err_mat_FL)
    disp(strcat('mean error:',num2str(mean(Outputs.Q.err_mat_FL))))
    ##
    Model.ShowComparison = copy('FullTimeSeries')
# ddd/Package_RunMe_Comparison.m:207
    Model.n_ts_samples = copy(100)
# ddd/Package_RunMe_Comparison.m:208
    Model.Guess_Q_TS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:209
    Model.Guess_Q_TS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_TS,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:210
    Model.Guess_Q_TS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_TS,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:211
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:212
    Model.Guess_Q_TS[arange(end() - Model.num_basis + 1,end())]=Model.Coeff_vecs[1]
# ddd/Package_RunMe_Comparison.m:213
    Outputs.Q.mean_err_TS,Grad_calc,Outputs.Q.err_mat_TS=Package_Q_Matrix_Comparisons(Model.Guess_Q_TS,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:214
    disp('=========================================================')
    disp('Full Time Series:')
    disp(Outputs.Q.err_mat_TS)
    disp(strcat('mean error:',num2str(mean(Outputs.Q.err_mat_TS))))
    ##
    Model.ShowComparison = copy('Composite')
# ddd/Package_RunMe_Comparison.m:221
    Model.Lambda = copy(0)
# ddd/Package_RunMe_Comparison.m:222
    Model.n_ts_samples = copy(10)
# ddd/Package_RunMe_Comparison.m:222
    Model.pts_of_int = copy(concat([1,4]))
# ddd/Package_RunMe_Comparison.m:222
    Model.Guess_Q_C,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:223
    ##
    Model.ShowComparison = copy('Composite')
# ddd/Package_RunMe_Comparison.m:225
    Model.Lambda = copy(0.01)
# ddd/Package_RunMe_Comparison.m:226
    Model.n_ts_samples = copy(10)
# ddd/Package_RunMe_Comparison.m:226
    Model.pts_of_int = copy(concat([1,4]))
# ddd/Package_RunMe_Comparison.m:226
    Model.Guess_Q_C,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_C,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:227
    ##
    Model.ShowComparison = copy('Composite')
# ddd/Package_RunMe_Comparison.m:229
    Model.Lambda = copy(0.01)
# ddd/Package_RunMe_Comparison.m:230
    Model.n_ts_samples = copy(10)
# ddd/Package_RunMe_Comparison.m:230
    Model.pts_of_int = copy(concat([1,4]))
# ddd/Package_RunMe_Comparison.m:230
    Model.Guess_Q_C,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_C,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:231
    ##
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:233
    #Model.Guess_Q_C(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
    Outputs.Q.mean_err_C,Grad_calc,Outputs.Q.err_mat_C=Package_Q_Matrix_Comparisons(Model.Guess_Q_C,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:235
    disp('=========================================================')
    disp('Composite:')
    disp(Outputs.Q.err_mat_C)
    disp(strcat('mean error:',num2str(mean(Outputs.Q.err_mat_C))))
    
    ##
    figure(1)
    plot(log10(1 + dot(89,Global.time_pts)),2 + log10(sqrt(Outputs.Q.err_mat_TS)),'--x','LineWidth',2)
    hold('on')
    plot(log10(1 + dot(89,Global.time_pts)),2 + log10(sqrt(Outputs.Q.err_mat_FS)),'--x','LineWidth',2)
    hold('on')
    plot(log10(1 + dot(89,Global.time_pts)),2 + log10(sqrt(Outputs.Q.err_mat_FL)),'--x','LineWidth',2)
    hold('on')
    plot(log10(1 + dot(89,Global.time_pts)),2 + log10(sqrt(Outputs.Q.err_mat_C)),'--x','LineWidth',2)
    hold('on')
    axis(concat([- 0.25,2,- 0.5,2.5]))
    legend('Trajectory','Full Snapshot','Two Snapshot','Integrated')
    #set(gca,'XTick',0:2)
#set(gca,'YTick',-3:3)
    set(gca,'fontsize',17)
    axis('square')
    tightfig()
    ##
    Model.Guess_P_FS = copy(Package_Func_QVec2PMat(Model.Guess_Q_FS,Model))
# ddd/Package_RunMe_Comparison.m:255
    Model.Guess_P_FL = copy(Package_Func_QVec2PMat(Model.Guess_Q_FL,Model))
# ddd/Package_RunMe_Comparison.m:256
    Model.Guess_P_TS = copy(Package_Func_QVec2PMat(Model.Guess_Q_TS,Model))
# ddd/Package_RunMe_Comparison.m:257
    Model.Guess_P_C = copy(Package_Func_QVec2PMat(Model.Guess_Q_C,Model))
# ddd/Package_RunMe_Comparison.m:258
    ##
    P_Lb=zeros(Model.num_basis)
# ddd/Package_RunMe_Comparison.m:261
    P_Lb[1 == eye(size(P_Lb))]=- Inf
# ddd/Package_RunMe_Comparison.m:262
    P_Lb[arange(),end() + 1]=0
# ddd/Package_RunMe_Comparison.m:263
    P_beq=sparse(Model.num_basis + 1,1)
# ddd/Package_RunMe_Comparison.m:264
    P_beq[end()]=1
# ddd/Package_RunMe_Comparison.m:265
    P_Aeq=sparse(Model.num_basis + 1,Model.num_basis ** 2)
# ddd/Package_RunMe_Comparison.m:266
    for kk in arange(1,(Model.num_basis + 1)).reshape(-1):
        P_Aeq[kk,arange(dot((kk - 1),Model.num_basis) + 1,dot(kk,Model.num_basis))]=1
# ddd/Package_RunMe_Comparison.m:268
    
    ##
    
    test_init_Q_vec=zeros(Model.num_basis,Model.num_basis + 1)
# ddd/Package_RunMe_Comparison.m:272
    test_init_Q_vec[arange(),end()]=Model.Coeff_vecs[1]
# ddd/Package_RunMe_Comparison.m:273
    ##
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:276
    Model.Guess_P_FS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:277
    Model.Guess_P_FS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_FS,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:278
    Model.Guess_P_FS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_FS,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:279
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:280
    #Model.Guess_P_FS(:,end) = Model.Coeff_vecs{1};
    Outputs.P.mean_err_FS,Grad_calc,Outputs.P.err_mat_FS=Package_P_Matrix_Comparisons(Model.Guess_P_FS,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:282
    disp('=========================================================')
    disp('Best we can hope for:')
    disp(Outputs.P.err_mat_FS)
    disp(strcat('mean error:',num2str(mean(Outputs.P.err_mat_FS))))
    ##
    Model.ShowComparison = copy('FirstLastPoint')
# ddd/Package_RunMe_Comparison.m:288
    Model.pts_of_int = copy(concat([1,4]))
# ddd/Package_RunMe_Comparison.m:289
    Model.Guess_P_FL,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:290
    Model.Guess_P_FL,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_FL,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:291
    Model.Guess_P_FL,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_FL,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:292
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:293
    #Model.Guess_P_FL(:,end) = Model.Coeff_vecs{1};
    Outputs.P.mean_err_FL,Grad_calc,Outputs.P.err_mat_FL=Package_P_Matrix_Comparisons(Model.Guess_P_FL,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:295
    disp('=========================================================')
    disp('Worst case scenario:')
    disp(Outputs.P.err_mat_FL)
    disp(strcat('mean error:',num2str(mean(Outputs.P.err_mat_FL))))
    ##
    Model.ShowComparison = copy('FullTimeSeries')
# ddd/Package_RunMe_Comparison.m:301
    Model.n_ts_samples = copy(100)
# ddd/Package_RunMe_Comparison.m:302
    Model.Guess_P_TS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:303
    Model.Guess_P_TS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_TS,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:304
    Model.Guess_P_TS,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_TS,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:305
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:306
    Model.Guess_P_TS[arange(),end()]=Model.Coeff_vecs[1]
# ddd/Package_RunMe_Comparison.m:307
    Outputs.P.mean_err_TS,Grad_calc,Outputs.P.err_mat_TS=Package_P_Matrix_Comparisons(Model.Guess_P_TS,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:308
    disp('=========================================================')
    disp('Full Time Series:')
    disp(Outputs.P.err_mat_TS)
    disp(strcat('mean error:',num2str(mean(Outputs.P.err_mat_TS))))
    ##
    Model.ShowComparison = copy('Composite')
# ddd/Package_RunMe_Comparison.m:315
    Model.Lambda = copy(0)
# ddd/Package_RunMe_Comparison.m:316
    Model.n_ts_samples = copy(10)
# ddd/Package_RunMe_Comparison.m:316
    Model.pts_of_int = copy(concat([1,4]))
# ddd/Package_RunMe_Comparison.m:316
    Model.Guess_P_C,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:317
    ##
    Model.ShowComparison = copy('Composite')
# ddd/Package_RunMe_Comparison.m:319
    Model.Lambda = copy(0.01)
# ddd/Package_RunMe_Comparison.m:320
    Model.n_ts_samples = copy(10)
# ddd/Package_RunMe_Comparison.m:320
    Model.pts_of_int = copy(concat([1,4]))
# ddd/Package_RunMe_Comparison.m:320
    Model.Guess_P_C,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_C,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:321
    ##
    Model.ShowComparison = copy('Composite')
# ddd/Package_RunMe_Comparison.m:323
    Model.Lambda = copy(0.01)
# ddd/Package_RunMe_Comparison.m:324
    Model.n_ts_samples = copy(10)
# ddd/Package_RunMe_Comparison.m:324
    Model.pts_of_int = copy(concat([1,4]))
# ddd/Package_RunMe_Comparison.m:324
    Model.Guess_P_C,__,tmp_exitflag,__,tmp_output=fmincon(lambda x=None: Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_C,[],[],P_Aeq,P_beq,P_Lb,[],[],opts,nargout=5)
# ddd/Package_RunMe_Comparison.m:325
    ##
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:327
    #Model.Guess_P_C(:,end) = Model.Coeff_vecs{1};
    Outputs.P.mean_err_C,Grad_calc,Outputs.P.err_mat_C=Package_P_Matrix_Comparisons(Model.Guess_P_C,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:329
    disp('=========================================================')
    disp('Composite:')
    disp(Outputs.P.err_mat_C)
    disp(strcat('mean error:',num2str(mean(Outputs.P.err_mat_C))))
    ##
    figure(2)
    plot(log10(1 + dot(89,Global.time_pts)),2 + log10(sqrt(Outputs.P.err_mat_TS)),'--x','LineWidth',2)
    hold('on')
    plot(log10(1 + dot(89,Global.time_pts)),2 + log10(sqrt(Outputs.P.err_mat_FS)),'--x','LineWidth',2)
    hold('on')
    plot(log10(1 + dot(89,Global.time_pts)),2 + log10(sqrt(Outputs.P.err_mat_FL)),'--x','LineWidth',2)
    hold('on')
    plot(log10(1 + dot(89,Global.time_pts)),2 + log10(sqrt(Outputs.P.err_mat_C)),'--x','LineWidth',2)
    hold('on')
    axis(concat([- 0.25,2,- 0.5,2.5]))
    legend('Trajectory','Full Snapshot','Two Snapshot','Integrated')
    #set(gca,'XTick',0:2)
#set(gca,'YTick',-3:3)
    set(gca,'fontsize',17)
    axis('square')
    tightfig()
    ##
    clc
    ##
    disp('All Time Series')
    disp(mean(sqrt(Outputs.Q.err_mat_TS(arange(2,end())))))
    disp(mean(sqrt(Outputs.P.err_mat_TS(arange(2,end())))))
    ##
    disp('All Snapshots')
    disp(mean(sqrt(Outputs.Q.err_mat_FS)))
    disp(mean(sqrt(Outputs.P.err_mat_FS)))
    ##
    disp('Two Snapshots')
    disp(mean(sqrt(Outputs.Q.err_mat_FL)))
    disp(mean(sqrt(Outputs.P.err_mat_FL)))
    ##
    disp('Composite')
    disp(mean(sqrt(Outputs.Q.err_mat_C)))
    disp(mean(sqrt(Outputs.P.err_mat_C)))
    ##
    
    ##
####### Checking Derivatives Q then P
    clc
    test_init_P_mat=Package_Func_QVec2PMat(test_init_Q_vec,Model)
# ddd/Package_RunMe_Comparison.m:378
    ##
    clc
    format('long')
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:382
    a,b,c=Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:383
    test_init_Q_vec[2]=test_init_Q_vec(2) + 0.0001
# ddd/Package_RunMe_Comparison.m:384
    A,B,C=Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:385
    disp(strcat('Derivative checking:',num2str(concat([(A - a) / 0.0001,b(2),B(2)]))))
    disp(strcat('Error:',num2str(a)))
    format('short')
    
    format('long')
    Model.ShowComparison = copy('FullSnapshot')
# ddd/Package_RunMe_Comparison.m:391
    a,b,__=Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:392
    test_init_P_mat[2]=test_init_P_mat(2) + 0.0001
# ddd/Package_RunMe_Comparison.m:393
    A,B,__=Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:394
    disp(strcat('Derivative checking:',num2str(concat([(A - a) / 0.0001,b(2),B(2)]))))
    disp(strcat('Error:',num2str(a)))
    format('short')
    ##
    clc
    format('long')
    Model.ShowComparison = copy('FirstLastPoint')
# ddd/Package_RunMe_Comparison.m:401
    Model.pts_of_int = copy(concat([1,6]))
# ddd/Package_RunMe_Comparison.m:402
    a,b,__=Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:403
    test_init_Q_vec[2]=test_init_Q_vec(2) + 0.0001
# ddd/Package_RunMe_Comparison.m:404
    A,B,__=Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:405
    disp(strcat('Derivative checking:',num2str(concat([(A - a) / 0.0001,b(2),B(2)]))))
    disp(strcat('Error:',num2str(a)))
    format('short')
    
    format('long')
    Model.ShowComparison = copy('FirstLastPoint')
# ddd/Package_RunMe_Comparison.m:411
    Model.pts_of_int = copy(concat([1,6]))
# ddd/Package_RunMe_Comparison.m:412
    a,b,__=Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:413
    test_init_P_mat[2]=test_init_P_mat(2) + 0.0001
# ddd/Package_RunMe_Comparison.m:414
    A,B,__=Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:415
    disp(strcat('Derivative checking:',num2str(concat([(A - a) / 0.0001,b(2),B(2)]))))
    disp(strcat('Error:',num2str(a)))
    format('short')
    ##
    clc
    format('long')
    Model.ShowComparison = copy('FullTimeSeries')
# ddd/Package_RunMe_Comparison.m:422
    Model.n_ts_samples = copy(100)
# ddd/Package_RunMe_Comparison.m:423
    a,b,__=Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:424
    test_init_Q_vec[2]=test_init_Q_vec(2) + 0.0001
# ddd/Package_RunMe_Comparison.m:425
    A,B,__=Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:426
    disp(strcat('Derivative checking:',num2str(concat([(A - a) / 0.0001,b(2),B(2)]))))
    disp(strcat('Error:',num2str(a)))
    format('short')
    
    format('long')
    Model.ShowComparison = copy('FullTimeSeries')
# ddd/Package_RunMe_Comparison.m:432
    Model.n_ts_samples = copy(100)
# ddd/Package_RunMe_Comparison.m:433
    a,b,__=Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:434
    test_init_P_mat[2]=test_init_P_mat(2) + 0.0001
# ddd/Package_RunMe_Comparison.m:435
    A,B,__=Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:436
    disp(strcat('Derivative checking:',num2str(concat([(A - a) / 0.0001,b(2),B(2)]))))
    disp(strcat('Error:',num2str(a)))
    format('short')
    ##
    clc
    format('long')
    Model.ShowComparison = copy('Composite')
# ddd/Package_RunMe_Comparison.m:443
    Model.Lambda = copy(0.25)
# ddd/Package_RunMe_Comparison.m:444
    Model.n_ts_samples = copy(10)
# ddd/Package_RunMe_Comparison.m:444
    Model.pts_of_int = copy(concat([1,6]))
# ddd/Package_RunMe_Comparison.m:444
    a,b,__=Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:445
    test_init_Q_vec[2]=test_init_Q_vec(2) + 0.0001
# ddd/Package_RunMe_Comparison.m:446
    A,B,__=Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:447
    disp(strcat('Derivative checking:',num2str(concat([(A - a) / 0.0001,b(2),B(2)]))))
    disp(strcat('Error:',num2str(a)))
    format('short')
    
    format('long')
    Model.ShowComparison = copy('Composite')
# ddd/Package_RunMe_Comparison.m:453
    Model.Lambda = copy(0.25)
# ddd/Package_RunMe_Comparison.m:454
    Model.n_ts_samples = copy(10)
# ddd/Package_RunMe_Comparison.m:454
    Model.pts_of_int = copy(concat([1,6]))
# ddd/Package_RunMe_Comparison.m:454
    a,b,__=Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:455
    test_init_P_mat[2]=test_init_P_mat(2) + 0.0001
# ddd/Package_RunMe_Comparison.m:456
    A,B,__=Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global,nargout=3)
# ddd/Package_RunMe_Comparison.m:457
    disp(strcat('Derivative checking:',num2str(concat([(A - a) / 0.0001,b(2),B(2)]))))
    disp(strcat('Error:',num2str(a)))
    format('short')