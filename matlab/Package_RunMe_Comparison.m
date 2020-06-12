clear all
close all
clc
% Load processed data
load('SDE_ProcessedData.mat')
% Rescale time to (0,1) interval
Global.time_pts = Global.time_pts./(Global.time_pts(end) - Global.time_pts(1));
Model.norm_fac = 'Relative';
%% Choose number and locations of basis function
Model.num_basis = 30;
[~,Model.Basis_Loc] = kmeans(Global.all_pts_MixModel,Model.num_basis,'Replicates',5);

%% Define form of basis functions used.
Model.subfunc.An = @(n,eps) 2*pi.^(n/2)./gamma(n/2);
Model.subfunc.Vn = @(n,eps) (eps.^n).*pi.^(n/2)./gamma((n/2)+1);
Model.subfunc.Intn = @(n,eps) eps.^n /(n*(n+1));
Model.basis_nd = @(x,n,eps) max( 0, 1-(vecnorm(x',2)./eps) )./(Model.subfunc.An(n,eps)*Model.subfunc.Intn(n,eps));

%% Choose scale parameter for each basis function. Connect to (minimum) N nearest basis functions.
Model.subfunc.connections = 2; % Minimum number of connections to other 
Model.Basis_eps = zeros(Model.num_basis, 1);
for ii = 1:Model.num_basis
    disp(ii)
    d_tmp = pdist2( Model.Basis_Loc(ii,:),Model.Basis_Loc );
    d_tmp = sort( d_tmp(d_tmp > 0));
    Model.Basis_eps(ii) = 1.00001*d_tmp(Model.subfunc.connections);
end
clear d_tmp ii 

%% Use Halton Node Placing for Monte Carlo integration
p = haltonset(Global.Data_Dim);
X0 = net(p,10^5);
for ii = 1:Model.num_basis
    disp(ii)
    for jj = ii:Model.num_basis
        tmp = Package_MC_Integration(ii, jj, Model, Global, X0);
        Model.MassMat(ii,jj) = tmp;
        Model.MassMat(jj,ii) = tmp;
    end
end
clear ii jj tmp X0 p 
Model.MassMat = sparse(Model.MassMat);

%% Plot Mass Matrix and Inverse
figure()
imagesc(log( Model.MassMat ))

Model.InvMassMat = inv(Model.MassMat);
Model.InvMassMat = (Model.InvMassMat + Model.InvMassMat')/2;

figure()
imagesc(log( abs( Model.InvMassMat ) ) )

%% Find coefficient vectors for each time pt.
Model.Coeff_vecs = cell(numel(Global.time_pts),1);
for vv = 1:numel(Global.time_pts)
    disp(vv)
    Model.Coeff_vecs{vv} = zeros(Model.num_basis,1);
    for ii = 1:Model.num_basis
    tmp = Model.basis_nd(Global.sample_pts{vv} - Model.Basis_Loc(ii,:),Global.Data_Dim, Model.Basis_eps(ii));
    Model.Coeff_vecs{vv}(ii) =  sum(tmp);
    end
    Model.Coeff_vecs{vv} = Model.Coeff_vecs{vv} / sum(Model.Coeff_vecs{vv});
end
clear vv ii tmp
%%
Model.TimeSeries = cell(numel(Global.time_pts),Global.Num_pts_each_time(1));
for ww = 1:Global.Num_pts_each_time(1)
    disp(ww)
    for vv = 1:numel(Global.time_pts)
        %disp(vv)
        Model.TimeSeries{vv,ww} = zeros(Model.num_basis,1);
        for ii = 1:Model.num_basis
            tmp = Model.basis_nd(Global.sample_pts{vv}(ww,:) - Model.Basis_Loc(ii,:),Global.Data_Dim, Model.Basis_eps(ii));
            Model.TimeSeries{vv,ww}(ii) =  tmp;
        end
        Model.TimeSeries{vv,ww} = Model.TimeSeries{vv,ww} / sum(Model.TimeSeries{vv,ww});
    end
end
clear vv ww ii tmp
%% Plot graphs of a few time points.
G_tmp = graph(Model.MassMat);
G_tmp = rmedge(G_tmp, 1:numnodes(G_tmp), 1:numnodes(G_tmp));
for vv = [1,9,11]
figure(vv),
plot(G_tmp,'XData',Model.Basis_Loc(:,1),'YData',Model.Basis_Loc(:,2)), hold on
scatter(Model.Basis_Loc(:,1),Model.Basis_Loc(:,2),50, Model.Coeff_vecs{vv},'filled')
axis([0 5 0 5])
end
clear G_tmp vv 
%% Plot graphs of a few time points.
G_tmp = graph(Model.MassMat);
G_tmp = rmedge(G_tmp, 1:numnodes(G_tmp), 1:numnodes(G_tmp));
for vv = 1:11
figure(vv),
plot(G_tmp,'XData',Model.Basis_Loc(:,1),'YData',Model.Basis_Loc(:,2)), hold on
scatter(Model.Basis_Loc(:,1),Model.Basis_Loc(:,2),50, Model.TimeSeries{vv,4},'filled')
axis([0 5 0 5])
end
clear G_tmp vv 


%% Plot coefficient matrices 
CoeffMat = cell2mat(Model.Coeff_vecs');
figure()
imagesc(CoeffMat )
drawnow

%%  Define useful quantities for later.
Model.InvSums = sum(Model.InvMassMat);
Model.Indc = sparse(Model.MassMat~=0);
Indc_sum = sum(Model.Indc);
Indc_zero = sparse(Model.MassMat==0);
Indc_vec = find(Model.MassMat~=0);
Indc_diag = find( ismember(  Indc_vec' , 1:( Model.num_basis+1 ): Model.num_basis^2) );


%% Define initial guess at P matrix
init_P = double(Model.Indc);
for vv = 1:Model.num_basis
    init_P(vv,vv) = 0;
    init_P(vv,vv) = - sum( init_P(:,vv) );
end
clear vv
%% Calculate corresponding Q and vectorise non-zero entries
init_Q = Model.MassMat*init_P;
init_Q(Indc_zero) = 0;
init_Q_vec = init_Q(~Indc_zero);
%% Specify effective dimension of problem.
eff_var = numel(init_Q_vec);
disp(strcat('Effective number of parameters:',num2str(eff_var)))
disp(strcat('Eff. basis functions for full DDD:',num2str(sqrt(eff_var))))
%% Append initial condition onto end of vector.
init_Q_vec(end+1:end+Model.num_basis) = Model.Coeff_vecs{1};
%% Creating constraints matrix for Aeq * init_Q_vec = beq
Q_Aeq = sparse(Model.num_basis + 1, Model.num_basis + eff_var);
run_tot_tmp = 1;
for kk = 1:Model.num_basis
    ind_tmp = run_tot_tmp:run_tot_tmp+Indc_sum(kk)-1;
    Q_Aeq(kk,ind_tmp) = Model.InvSums(Model.Indc(kk,:));
    run_tot_tmp = run_tot_tmp + Indc_sum(kk);
end
clear kk ind_tmp run_tot_tmp
Q_Aeq( Model.num_basis + 1 , ( end-Model.num_basis+1 ):end) = 1;
Q_beq = sparse(Model.num_basis + 1,1);
Q_beq(end) = 1;
%% 
imagesc ( [[sum(Model.InvMassMat*init_Q)';nan], Q_Aeq * init_Q_vec, Q_beq] )

%% Create positivity constraint Alt * init_Q_vec <= blt
Q_Alt = sparse(Model.num_basis + eff_var, Model.num_basis + eff_var);
run_tot_tmp = 1;
for kk = 1:Model.num_basis
    xx = find( Model.Indc(kk,:) );
    yy = Model.InvMassMat(xx,xx);
    ind_tmp = run_tot_tmp:run_tot_tmp+Indc_sum(kk)-1;
    Q_Alt(ind_tmp,ind_tmp)  = -yy; 
    run_tot_tmp = run_tot_tmp + Indc_sum(kk);
end
ind_tmp = run_tot_tmp:run_tot_tmp+Model.num_basis-1;
Q_Alt(ind_tmp,ind_tmp) = -eye(Model.num_basis);
for kk = Indc_diag
    Q_Alt(kk,:) = - Q_Alt(kk,:);
end
Q_blt = sparse(Model.num_basis + eff_var,1);
clear kk ind_tmp
%%

opts = optimoptions('fmincon','Display','iter-detailed');
opts.UseParallel = true(1,1);
opts.MaxIterations = 3000;
opts.SpecifyObjectiveGradient = true(1,1);
opts.Algorithm = 'sqp';
%opts.TypicalX = init_K;
opts.MaxIterations = 300;
%%
%[Model.Guess_Q,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Find_Q_Matrix(x,Model,Global),init_Q_vec,Alt,blt,Aeq,beq,[],[],[],opts);
%%
test_init_Q_vec = init_Q_vec;
test_init_Q_vec(1:eff_var) = 0;
%%
Model.ShowComparison = 'FullSnapshot';
[Model.Guess_Q_FS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
[Model.Guess_Q_FS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_FS,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
[Model.Guess_Q_FS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_FS,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
Model.ShowComparison = 'FullSnapshot';
%Model.Guess_Q_FS(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
[Outputs.Q.mean_err_FS,Grad_calc,Outputs.Q.err_mat_FS] = Package_Q_Matrix_Comparisons(Model.Guess_Q_FS,Model,Global);
disp('=========================================================')
disp('Best we can hope for:')
disp(Outputs.Q.err_mat_FS)
disp(strcat('mean error:',num2str(mean(Outputs.Q.err_mat_FS))))
%%
Model.ShowComparison = 'FirstLastPoint';
Model.pts_of_int = [1,4];
[Model.Guess_Q_FL,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
[Model.Guess_Q_FL,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_FL,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
[Model.Guess_Q_FL,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_FL,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
Model.ShowComparison = 'FullSnapshot';
%Model.Guess_Q_FL(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
[Outputs.Q.mean_err_FL,Grad_calc,Outputs.Q.err_mat_FL] = Package_Q_Matrix_Comparisons(Model.Guess_Q_FL,Model,Global);
disp('=========================================================')
disp('Worst case scenario:')
disp(Outputs.Q.err_mat_FL)
disp(strcat('mean error:',num2str(mean(Outputs.Q.err_mat_FL))))
%%
Model.ShowComparison = 'FullTimeSeries';
Model.n_ts_samples = 100;
[Model.Guess_Q_TS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
[Model.Guess_Q_TS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_TS,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
[Model.Guess_Q_TS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_TS,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
Model.ShowComparison = 'FullSnapshot';
Model.Guess_Q_TS(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
[Outputs.Q.mean_err_TS,Grad_calc,Outputs.Q.err_mat_TS] = Package_Q_Matrix_Comparisons(Model.Guess_Q_TS,Model,Global);
disp('=========================================================')
disp('Full Time Series:')
disp(Outputs.Q.err_mat_TS)
disp(strcat('mean error:',num2str(mean(Outputs.Q.err_mat_TS))))

%%
Model.ShowComparison = 'Composite';
Model.Lambda = 0; Model.n_ts_samples = 10; Model.pts_of_int = [1,4];
[Model.Guess_Q_C,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
%%
Model.ShowComparison = 'Composite';
Model.Lambda = 0.01; Model.n_ts_samples = 10; Model.pts_of_int = [1,4];
[Model.Guess_Q_C,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_C,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
%%
Model.ShowComparison = 'Composite';
Model.Lambda = 0.01; Model.n_ts_samples = 10; Model.pts_of_int = [1,4];
[Model.Guess_Q_C,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_Q_Matrix_Comparisons(x,Model,Global),Model.Guess_Q_C,Q_Alt,Q_blt,Q_Aeq,Q_beq,[],[],[],opts);
%%
Model.ShowComparison = 'FullSnapshot';
%Model.Guess_Q_C(end-Model.num_basis+1:end) = Model.Coeff_vecs{1};
[Outputs.Q.mean_err_C,Grad_calc,Outputs.Q.err_mat_C] = Package_Q_Matrix_Comparisons(Model.Guess_Q_C,Model,Global);
disp('=========================================================')
disp('Composite:')
disp(Outputs.Q.err_mat_C)
disp(strcat('mean error:',num2str(mean(Outputs.Q.err_mat_C))))
 
%%
figure(1)
plot( log10( 1 + 89*Global.time_pts) , 2 + log10( sqrt( Outputs.Q.err_mat_TS )) , '--x','LineWidth',2 ), hold on
plot( log10( 1 + 89*Global.time_pts) , 2 + log10( sqrt( Outputs.Q.err_mat_FS )) , '--x','LineWidth',2 ), hold on
plot( log10( 1 + 89*Global.time_pts) , 2 + log10( sqrt( Outputs.Q.err_mat_FL )) , '--x','LineWidth',2  ), hold on
plot( log10( 1 + 89*Global.time_pts) , 2 + log10( sqrt( Outputs.Q.err_mat_C )) , '--x','LineWidth',2  ), hold on
axis([-0.25, 2, -0.5, 2.5])
legend('Trajectory', 'Full Snapshot', 'Two Snapshot', 'Integrated')
%set(gca,'XTick',0:2)
%set(gca,'YTick',-3:3)
   set(gca,'fontsize',17)
axis square
 tightfig(); 
%%
Model.Guess_P_FS = Package_Func_QVec2PMat(Model.Guess_Q_FS, Model);
Model.Guess_P_FL = Package_Func_QVec2PMat(Model.Guess_Q_FL, Model);
Model.Guess_P_TS = Package_Func_QVec2PMat(Model.Guess_Q_TS, Model);
Model.Guess_P_C = Package_Func_QVec2PMat(Model.Guess_Q_C, Model);

%%
P_Lb = zeros(Model.num_basis);
P_Lb( 1 == eye(size(P_Lb)) ) = -Inf;
P_Lb(:,end+1) = 0;
P_beq = sparse(Model.num_basis + 1,1);
P_beq(end) = 1;
P_Aeq = sparse(Model.num_basis + 1, Model.num_basis^2);
for kk = 1:(Model.num_basis + 1)
    P_Aeq(kk,(kk-1)*Model.num_basis+1:kk*Model.num_basis) = 1;
end
%%

test_init_Q_vec = zeros(Model.num_basis, Model.num_basis + 1);
test_init_Q_vec(:,end) = Model.Coeff_vecs{1};

%%
Model.ShowComparison = 'FullSnapshot';
[Model.Guess_P_FS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
[Model.Guess_P_FS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_FS,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
[Model.Guess_P_FS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_FS,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
Model.ShowComparison = 'FullSnapshot';
%Model.Guess_P_FS(:,end) = Model.Coeff_vecs{1};
[Outputs.P.mean_err_FS,Grad_calc,Outputs.P.err_mat_FS] = Package_P_Matrix_Comparisons(Model.Guess_P_FS,Model,Global);
disp('=========================================================')
disp('Best we can hope for:')
disp(Outputs.P.err_mat_FS)
disp(strcat('mean error:',num2str(mean(Outputs.P.err_mat_FS))))
%%
Model.ShowComparison = 'FirstLastPoint';
Model.pts_of_int = [1,4];
[Model.Guess_P_FL,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
[Model.Guess_P_FL,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_FL,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
[Model.Guess_P_FL,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_FL,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
Model.ShowComparison = 'FullSnapshot';
%Model.Guess_P_FL(:,end) = Model.Coeff_vecs{1};
[Outputs.P.mean_err_FL,Grad_calc,Outputs.P.err_mat_FL] = Package_P_Matrix_Comparisons(Model.Guess_P_FL,Model,Global);
disp('=========================================================')
disp('Worst case scenario:')
disp(Outputs.P.err_mat_FL)
disp(strcat('mean error:',num2str(mean(Outputs.P.err_mat_FL))))
%%
Model.ShowComparison = 'FullTimeSeries';
Model.n_ts_samples = 100;
[Model.Guess_P_TS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
[Model.Guess_P_TS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_TS,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
[Model.Guess_P_TS,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_TS,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
Model.ShowComparison = 'FullSnapshot';
Model.Guess_P_TS(:,end) = Model.Coeff_vecs{1};
[Outputs.P.mean_err_TS,Grad_calc,Outputs.P.err_mat_TS] = Package_P_Matrix_Comparisons(Model.Guess_P_TS,Model,Global);
disp('=========================================================')
disp('Full Time Series:')
disp(Outputs.P.err_mat_TS)
disp(strcat('mean error:',num2str(mean(Outputs.P.err_mat_TS))))

%%
Model.ShowComparison = 'Composite';
Model.Lambda = 0; Model.n_ts_samples = 10; Model.pts_of_int = [1,4];
[Model.Guess_P_C,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),test_init_Q_vec,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
%%
Model.ShowComparison = 'Composite';
Model.Lambda = 0.01; Model.n_ts_samples = 10; Model.pts_of_int = [1,4];
[Model.Guess_P_C,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_C,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
%%
Model.ShowComparison = 'Composite';
Model.Lambda = 0.01; Model.n_ts_samples = 10; Model.pts_of_int = [1,4];
[Model.Guess_P_C,~,tmp_exitflag,~,tmp_output] = fmincon(@(x)Package_P_Matrix_Comparisons(x,Model,Global),Model.Guess_P_C,[],[],P_Aeq,P_beq,P_Lb,[],[],opts);
%%
Model.ShowComparison = 'FullSnapshot';
%Model.Guess_P_C(:,end) = Model.Coeff_vecs{1};
[Outputs.P.mean_err_C,Grad_calc,Outputs.P.err_mat_C] = Package_P_Matrix_Comparisons(Model.Guess_P_C,Model,Global);
disp('=========================================================')
disp('Composite:')
disp(Outputs.P.err_mat_C)
disp(strcat('mean error:',num2str(mean(Outputs.P.err_mat_C))))


%%
figure(2)
plot( log10( 1 + 89*Global.time_pts) , 2 + log10( sqrt(Outputs.P.err_mat_TS) ) , '--x','LineWidth',2 ), hold on
plot( log10( 1 + 89*Global.time_pts) , 2 + log10( sqrt(Outputs.P.err_mat_FS) ) , '--x','LineWidth',2 ), hold on
plot( log10( 1 + 89*Global.time_pts) , 2 + log10( sqrt(Outputs.P.err_mat_FL) ) , '--x','LineWidth',2  ), hold on
plot( log10( 1 + 89*Global.time_pts) , 2 + log10( sqrt(Outputs.P.err_mat_C) ) , '--x','LineWidth',2  ), hold on
axis([-0.25, 2, -0.5, 2.5])
legend('Trajectory', 'Full Snapshot', 'Two Snapshot', 'Integrated')
%set(gca,'XTick',0:2)
%set(gca,'YTick',-3:3)
   set(gca,'fontsize',17)
axis square
tightfig();

%%
clc
%%
disp("All Time Series")
disp( mean ( sqrt( Outputs.Q.err_mat_TS(2:end) ) ) )
disp( mean ( sqrt( Outputs.P.err_mat_TS(2:end) ) ) )
%%
disp("All Snapshots")
disp( mean ( sqrt( Outputs.Q.err_mat_FS ) ) )
disp( mean ( sqrt( Outputs.P.err_mat_FS ) ) )
%%
disp("Two Snapshots")
disp( mean ( sqrt( Outputs.Q.err_mat_FL ) ) )
disp( mean ( sqrt( Outputs.P.err_mat_FL ) ) )
%%
disp("Composite")
disp( mean ( sqrt( Outputs.Q.err_mat_C ) ) )
disp( mean ( sqrt( Outputs.P.err_mat_C ) ) )
%%






%%
%%%%%%% Checking Derivatives Q then P
clc
test_init_P_mat = Package_Func_QVec2PMat(test_init_Q_vec, Model);
%%
clc
format long
Model.ShowComparison = 'FullSnapshot';
[a,b,c] = Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global);
test_init_Q_vec(2) = test_init_Q_vec(2) + 1e-4;
[A,B,C] = Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global);
disp(strcat("Derivative checking:",num2str( [ (A - a)/1e-4, b(2), B(2) ] )))
disp(strcat("Error:",num2str( a )))
format short
%
format long
Model.ShowComparison = 'FullSnapshot';
[a,b,~] = Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global);
test_init_P_mat(2) = test_init_P_mat(2) + 1e-4;
[A,B,~] = Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global);
disp(strcat("Derivative checking:",num2str( [ (A - a)/1e-4, b(2), B(2) ] )))
disp(strcat("Error:",num2str( a )))
format short
%%
clc
format long
Model.ShowComparison = 'FirstLastPoint';
Model.pts_of_int = [1,6];
[a,b,~] = Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global);
test_init_Q_vec(2) = test_init_Q_vec(2) + 1e-4;
[A,B,~] = Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global);
disp(strcat("Derivative checking:",num2str( [ (A - a)/1e-4, b(2), B(2) ] )))
disp(strcat("Error:",num2str( a )))
format short
%
format long
Model.ShowComparison = 'FirstLastPoint';
Model.pts_of_int = [1,6];
[a,b,~] = Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global);
test_init_P_mat(2) = test_init_P_mat(2) + 1e-4;
[A,B,~] = Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global);
disp(strcat("Derivative checking:",num2str( [ (A - a)/1e-4, b(2), B(2) ] )))
disp(strcat("Error:",num2str( a )))
format short
%%
clc
format long
Model.ShowComparison = 'FullTimeSeries';
Model.n_ts_samples = 100;
[a,b,~] = Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global);
test_init_Q_vec(2) = test_init_Q_vec(2) + 1e-4;
[A,B,~] = Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global);
disp(strcat("Derivative checking:",num2str( [ (A - a)/1e-4, b(2), B(2) ] )))
disp(strcat("Error:",num2str( a )))
format short
%
format long
Model.ShowComparison = 'FullTimeSeries';
Model.n_ts_samples = 100;
[a,b,~] = Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global);
test_init_P_mat(2) = test_init_P_mat(2) + 1e-4;
[A,B,~] = Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global);
disp(strcat("Derivative checking:",num2str( [ (A - a)/1e-4, b(2), B(2) ] )))
disp(strcat("Error:",num2str( a )))
format short
%%
clc
format long
Model.ShowComparison = 'Composite';
Model.Lambda = 0.25; Model.n_ts_samples = 10; Model.pts_of_int = [1,6];
[a,b,~] = Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global);
test_init_Q_vec(2) = test_init_Q_vec(2) + 1e-4;
[A,B,~] = Package_Q_Matrix_Comparisons(test_init_Q_vec,Model,Global);
disp(strcat("Derivative checking:",num2str( [ (A - a)/1e-4, b(2), B(2) ] )))
disp(strcat("Error:",num2str( a )))
format short
%
format long
Model.ShowComparison = 'Composite';
Model.Lambda = 0.25; Model.n_ts_samples = 10; Model.pts_of_int = [1,6];
[a,b,~] = Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global);
test_init_P_mat(2) = test_init_P_mat(2) + 1e-4;
[A,B,~] = Package_P_Matrix_Comparisons(test_init_P_mat,Model,Global);
disp(strcat("Derivative checking:",num2str( [ (A - a)/1e-4, b(2), B(2) ] )))
disp(strcat("Error:",num2str( a )))
format short