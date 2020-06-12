function [error,Grad_mat,err_mat] = Func_P_Matrix_Comparisons(P_Mat,Model,Global)
tic

Start_Vec = P_Mat(:,end);
PF_Mat = P_Mat(:,1:end-1);



%% Full time
if Model.ShowComparison == "FullSnapshot"
    err_mat = zeros(numel(Global.time_pts),1);
    Grad_mat = zeros(size(PF_Mat,1),size(PF_Mat,1)+1);
    for vv = 1:numel(Global.time_pts)
        t_gap = (Global.time_pts(vv)-Global.time_pts(1));
        Coeff_vec = Model.Coeff_vecs{vv};
        [error_tmp,Grad_mat_tmp] = Func_SingleTimePointError(PF_Mat,Coeff_vec,Start_Vec,t_gap,Model);
        err_mat(vv) = error_tmp;
        Grad_mat = Grad_mat + Grad_mat_tmp;
    end
    
    error = sum(err_mat)/Global.num_tpts;
    Grad_mat = Grad_mat/Global.num_tpts;
    
    
    
    %disp(error)
end

%% Evaluate first and last point only!
if Model.ShowComparison == "FirstLastPoint"
    err_mat = zeros(2,1);
    Grad_mat = zeros(size(PF_Mat,1),size(PF_Mat,1)+1);
    for vv = Model.pts_of_int
        t_gap = (Global.time_pts(vv)-Global.time_pts(1));
        Coeff_vec = Model.Coeff_vecs{vv};
        [error_tmp,Grad_mat_tmp] = Func_SingleTimePointError(PF_Mat,Coeff_vec,Start_Vec,t_gap,Model);
        err_mat(vv) = error_tmp;
        Grad_mat = Grad_mat + Grad_mat_tmp;
    end
    err_mat = err_mat(Model.pts_of_int);
    error = sum(err_mat)/2;
    Grad_mat = Grad_mat/2;
    
    
    %disp(error)
end


%% Time Series
if Model.ShowComparison == "FullTimeSeries"
    
    err_mat = zeros( [size( Model.TimeSeries,1)-1,Model.n_ts_samples] );
    Grad_mat = zeros(size(PF_Mat,1),size(PF_Mat,1)+1);
    tot_valid = 0;
    for ww = 1:Model.n_ts_samples
    %for ww = randsample(1000,Model.n_ts_samples)'
        err_vec = zeros(1,numel(Global.time_pts)-1);
        Grad_mat_cache = zeros(size(PF_Mat,1),size(PF_Mat,1)+1);
        for vv = 2:numel(Global.time_pts)
            t_gap = (Global.time_pts(vv)-Global.time_pts(1));
            Coeff_vec = Model.TimeSeries{vv,ww};
            [error_tmp,Grad_mat_tmp] = Func_SingleTimePointError(PF_Mat,Coeff_vec,Model.TimeSeries{1,ww},t_gap,Model);
            err_vec(vv-1) = error_tmp;
            Grad_mat_cache = Grad_mat_cache + Grad_mat_tmp;
        end
        if sum( isnan(err_vec) ) == 0
            err_mat(:,ww) = err_vec;
            Grad_mat = Grad_mat + Grad_mat_cache;
            tot_valid = tot_valid + 1;
        end
    end
    
    err_mat(:,sum(err_mat) == 0) = [];
    
    error = sum(sum(err_mat))/((Global.num_tpts-1)*tot_valid);
    Grad_mat = Grad_mat/((Global.num_tpts-1)*tot_valid);
    
    %disp(error)
    
    %Grad_mat(:,1:Model.num_basis) = Model.MassMat\Grad_mat(:,1:Model.num_basis);
    Grad_mat(:,end) = 0;
    
    %Grad_mat = Grad_mat(logical( [Model.Indc,ones(Model.num_basis,1)]) );
    
end

%% Comparison
if Model.ShowComparison == "Composite"
    
    if Model.Lambda ~= 0
        
        err_mat_TS = zeros( [size( Model.TimeSeries,1)-1,Model.n_ts_samples] );
        Grad_mat_TS = zeros(size(PF_Mat,1),size(PF_Mat,1)+1);
        tot_valid = 0;
        for ww = 1:Model.n_ts_samples
            err_vec = zeros(1,numel(Global.time_pts)-1);
            Grad_mat_cache = zeros(size(PF_Mat,1),size(PF_Mat,1)+1);
            for vv = 2:numel(Global.time_pts)
                t_gap = (Global.time_pts(vv)-Global.time_pts(1));
                Coeff_vec = Model.TimeSeries{vv,ww};
                [error_tmp,Grad_mat_tmp] = Func_SingleTimePointError(PF_Mat,Coeff_vec,Model.TimeSeries{1,ww},t_gap,Model);
                err_vec(vv-1) = error_tmp;
                Grad_mat_cache = Grad_mat_cache + Grad_mat_tmp;
            end
            if sum( isnan(err_vec) ) == 0
                err_mat_TS(:,ww) = err_vec;
                Grad_mat_TS = Grad_mat_TS + Grad_mat_cache;
                tot_valid = tot_valid + 1;
            end
        end
        err_mat_TS(:,sum(err_mat_TS) == 0) = [];
        
        error_TS = sum(sum(err_mat_TS))/((Global.num_tpts-1)*tot_valid);
        Grad_mat_TS = Grad_mat_TS/((Global.num_tpts-1)*tot_valid);
        
        Grad_mat_TS(:,end) = 0;
        
    end
    
    % First last Snapshot
    
    err_mat_FL = zeros(2,1);
    Grad_mat_FL = zeros(size(PF_Mat,1),size(PF_Mat,1)+1);
    for vv = Model.pts_of_int
        t_gap = (Global.time_pts(vv)-Global.time_pts(1));
        Coeff_vec = Model.Coeff_vecs{vv};
        [error_tmp,Grad_mat_tmp] = Func_SingleTimePointError(PF_Mat,Coeff_vec,Start_Vec,t_gap,Model);
        err_mat_FL(vv) = error_tmp;
        Grad_mat_FL = Grad_mat_FL + Grad_mat_tmp;
    end
    err_mat_FL = err_mat_FL(Model.pts_of_int);
    error_FL = sum(err_mat_FL)/2;
    Grad_mat_FL = Grad_mat_FL/2;
    
    % Add up terms
    if Model.Lambda ~= 0
        error = Model.Lambda*error_TS + (1 - Model.Lambda)*error_FL;
        Grad_mat = Model.Lambda*Grad_mat_TS + (1 - Model.Lambda)*Grad_mat_FL;
        disp([error_TS, error_FL, error])
    else
        error = error_FL;
        Grad_mat = Grad_mat_FL;
        %disp(error)
    end
    
    err_mat = [];
    
    % Vectorise and remove bad entries
    
    
    
    
end

Grad_mat(abs( Grad_mat) > 10^5 ) = 0;

%disp(error)

%% time keeping
stp_time = toc;
if stp_time > 5
    disp(strcat('slow eval:',num2str(stp_time)))
end

end
