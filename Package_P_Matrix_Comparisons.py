# Generated with SMOP  0.41-beta
from libsmop import *
# ddd/Package_P_Matrix_Comparisons.m

    
@function
def Package_P_Matrix_Comparisons(P_Mat=None,Model=None,Global=None,*args,**kwargs):
    varargin = Package_P_Matrix_Comparisons.varargin
    nargin = Package_P_Matrix_Comparisons.nargin

    tic
    Start_Vec=P_Mat(arange(),end())
# ddd/Package_P_Matrix_Comparisons.m:4
    PF_Mat=P_Mat(arange(),arange(1,end() - 1))
# ddd/Package_P_Matrix_Comparisons.m:5
    ## Full time
    if Model.ShowComparison == 'FullSnapshot':
        err_mat=zeros(numel(Global.time_pts),1)
# ddd/Package_P_Matrix_Comparisons.m:11
        Grad_mat=zeros(size(PF_Mat,1),size(PF_Mat,1) + 1)
# ddd/Package_P_Matrix_Comparisons.m:12
        for vv in arange(1,numel(Global.time_pts)).reshape(-1):
            t_gap=(Global.time_pts(vv) - Global.time_pts(1))
# ddd/Package_P_Matrix_Comparisons.m:14
            Coeff_vec=Model.Coeff_vecs[vv]
# ddd/Package_P_Matrix_Comparisons.m:15
            error_tmp,Grad_mat_tmp=Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Start_Vec,t_gap,Model,nargout=2)
# ddd/Package_P_Matrix_Comparisons.m:16
            err_mat[vv]=error_tmp
# ddd/Package_P_Matrix_Comparisons.m:17
            Grad_mat=Grad_mat + Grad_mat_tmp
# ddd/Package_P_Matrix_Comparisons.m:18
        error=sum(err_mat) / Global.num_tpts
# ddd/Package_P_Matrix_Comparisons.m:21
        Grad_mat=Grad_mat / Global.num_tpts
# ddd/Package_P_Matrix_Comparisons.m:22
    
    ## Evaluate first and last point only!
    if Model.ShowComparison == 'FirstLastPoint':
        err_mat=zeros(2,1)
# ddd/Package_P_Matrix_Comparisons.m:31
        Grad_mat=zeros(size(PF_Mat,1),size(PF_Mat,1) + 1)
# ddd/Package_P_Matrix_Comparisons.m:32
        for vv in Model.pts_of_int.reshape(-1):
            t_gap=(Global.time_pts(vv) - Global.time_pts(1))
# ddd/Package_P_Matrix_Comparisons.m:34
            Coeff_vec=Model.Coeff_vecs[vv]
# ddd/Package_P_Matrix_Comparisons.m:35
            error_tmp,Grad_mat_tmp=Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Start_Vec,t_gap,Model,nargout=2)
# ddd/Package_P_Matrix_Comparisons.m:36
            err_mat[vv]=error_tmp
# ddd/Package_P_Matrix_Comparisons.m:37
            Grad_mat=Grad_mat + Grad_mat_tmp
# ddd/Package_P_Matrix_Comparisons.m:38
        err_mat=err_mat(Model.pts_of_int)
# ddd/Package_P_Matrix_Comparisons.m:40
        error=sum(err_mat) / 2
# ddd/Package_P_Matrix_Comparisons.m:41
        Grad_mat=Grad_mat / 2
# ddd/Package_P_Matrix_Comparisons.m:42
    
    ## Time Series
    if Model.ShowComparison == 'FullTimeSeries':
        err_mat=zeros(concat([size(Model.TimeSeries,1) - 1,Model.n_ts_samples]))
# ddd/Package_P_Matrix_Comparisons.m:52
        Grad_mat=zeros(size(PF_Mat,1),size(PF_Mat,1) + 1)
# ddd/Package_P_Matrix_Comparisons.m:53
        tot_valid=0
# ddd/Package_P_Matrix_Comparisons.m:54
        for ww in arange(1,Model.n_ts_samples).reshape(-1):
            #for ww = randsample(1000,Model.n_ts_samples)'
            err_vec=zeros(1,numel(Global.time_pts) - 1)
# ddd/Package_P_Matrix_Comparisons.m:57
            Grad_mat_cache=zeros(size(PF_Mat,1),size(PF_Mat,1) + 1)
# ddd/Package_P_Matrix_Comparisons.m:58
            for vv in arange(2,numel(Global.time_pts)).reshape(-1):
                t_gap=(Global.time_pts(vv) - Global.time_pts(1))
# ddd/Package_P_Matrix_Comparisons.m:60
                Coeff_vec=Model.TimeSeries[vv,ww]
# ddd/Package_P_Matrix_Comparisons.m:61
                error_tmp,Grad_mat_tmp=Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Model.TimeSeries[1,ww],t_gap,Model,nargout=2)
# ddd/Package_P_Matrix_Comparisons.m:62
                err_vec[vv - 1]=error_tmp
# ddd/Package_P_Matrix_Comparisons.m:63
                Grad_mat_cache=Grad_mat_cache + Grad_mat_tmp
# ddd/Package_P_Matrix_Comparisons.m:64
            if sum(isnan(err_vec)) == 0:
                err_mat[arange(),ww]=err_vec
# ddd/Package_P_Matrix_Comparisons.m:67
                Grad_mat=Grad_mat + Grad_mat_cache
# ddd/Package_P_Matrix_Comparisons.m:68
                tot_valid=tot_valid + 1
# ddd/Package_P_Matrix_Comparisons.m:69
        err_mat[arange(),sum(err_mat) == 0]=[]
# ddd/Package_P_Matrix_Comparisons.m:73
        error=sum(sum(err_mat)) / (dot((Global.num_tpts - 1),tot_valid))
# ddd/Package_P_Matrix_Comparisons.m:75
        Grad_mat=Grad_mat / (dot((Global.num_tpts - 1),tot_valid))
# ddd/Package_P_Matrix_Comparisons.m:76
        #Grad_mat(:,1:Model.num_basis) = Model.MassMat\Grad_mat(:,1:Model.num_basis);
        Grad_mat[arange(),end()]=0
# ddd/Package_P_Matrix_Comparisons.m:81
    
    ## Comparison
    if Model.ShowComparison == 'Composite':
        if Model.Lambda != 0:
            err_mat_TS=zeros(concat([size(Model.TimeSeries,1) - 1,Model.n_ts_samples]))
# ddd/Package_P_Matrix_Comparisons.m:92
            Grad_mat_TS=zeros(size(PF_Mat,1),size(PF_Mat,1) + 1)
# ddd/Package_P_Matrix_Comparisons.m:93
            tot_valid=0
# ddd/Package_P_Matrix_Comparisons.m:94
            for ww in arange(1,Model.n_ts_samples).reshape(-1):
                err_vec=zeros(1,numel(Global.time_pts) - 1)
# ddd/Package_P_Matrix_Comparisons.m:96
                Grad_mat_cache=zeros(size(PF_Mat,1),size(PF_Mat,1) + 1)
# ddd/Package_P_Matrix_Comparisons.m:97
                for vv in arange(2,numel(Global.time_pts)).reshape(-1):
                    t_gap=(Global.time_pts(vv) - Global.time_pts(1))
# ddd/Package_P_Matrix_Comparisons.m:99
                    Coeff_vec=Model.TimeSeries[vv,ww]
# ddd/Package_P_Matrix_Comparisons.m:100
                    error_tmp,Grad_mat_tmp=Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Model.TimeSeries[1,ww],t_gap,Model,nargout=2)
# ddd/Package_P_Matrix_Comparisons.m:101
                    err_vec[vv - 1]=error_tmp
# ddd/Package_P_Matrix_Comparisons.m:102
                    Grad_mat_cache=Grad_mat_cache + Grad_mat_tmp
# ddd/Package_P_Matrix_Comparisons.m:103
                if sum(isnan(err_vec)) == 0:
                    err_mat_TS[arange(),ww]=err_vec
# ddd/Package_P_Matrix_Comparisons.m:106
                    Grad_mat_TS=Grad_mat_TS + Grad_mat_cache
# ddd/Package_P_Matrix_Comparisons.m:107
                    tot_valid=tot_valid + 1
# ddd/Package_P_Matrix_Comparisons.m:108
            err_mat_TS[arange(),sum(err_mat_TS) == 0]=[]
# ddd/Package_P_Matrix_Comparisons.m:111
            error_TS=sum(sum(err_mat_TS)) / (dot((Global.num_tpts - 1),tot_valid))
# ddd/Package_P_Matrix_Comparisons.m:113
            Grad_mat_TS=Grad_mat_TS / (dot((Global.num_tpts - 1),tot_valid))
# ddd/Package_P_Matrix_Comparisons.m:114
            Grad_mat_TS[arange(),end()]=0
# ddd/Package_P_Matrix_Comparisons.m:116
        # First last Snapshot
        err_mat_FL=zeros(2,1)
# ddd/Package_P_Matrix_Comparisons.m:122
        Grad_mat_FL=zeros(size(PF_Mat,1),size(PF_Mat,1) + 1)
# ddd/Package_P_Matrix_Comparisons.m:123
        for vv in Model.pts_of_int.reshape(-1):
            t_gap=(Global.time_pts(vv) - Global.time_pts(1))
# ddd/Package_P_Matrix_Comparisons.m:125
            Coeff_vec=Model.Coeff_vecs[vv]
# ddd/Package_P_Matrix_Comparisons.m:126
            error_tmp,Grad_mat_tmp=Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Start_Vec,t_gap,Model,nargout=2)
# ddd/Package_P_Matrix_Comparisons.m:127
            err_mat_FL[vv]=error_tmp
# ddd/Package_P_Matrix_Comparisons.m:128
            Grad_mat_FL=Grad_mat_FL + Grad_mat_tmp
# ddd/Package_P_Matrix_Comparisons.m:129
        err_mat_FL=err_mat_FL(Model.pts_of_int)
# ddd/Package_P_Matrix_Comparisons.m:131
        error_FL=sum(err_mat_FL) / 2
# ddd/Package_P_Matrix_Comparisons.m:132
        Grad_mat_FL=Grad_mat_FL / 2
# ddd/Package_P_Matrix_Comparisons.m:133
        if Model.Lambda != 0:
            error=dot(Model.Lambda,error_TS) + dot((1 - Model.Lambda),error_FL)
# ddd/Package_P_Matrix_Comparisons.m:137
            Grad_mat=dot(Model.Lambda,Grad_mat_TS) + dot((1 - Model.Lambda),Grad_mat_FL)
# ddd/Package_P_Matrix_Comparisons.m:138
            disp(concat([error_TS,error_FL,error]))
        else:
            error=copy(error_FL)
# ddd/Package_P_Matrix_Comparisons.m:141
            Grad_mat=copy(Grad_mat_FL)
# ddd/Package_P_Matrix_Comparisons.m:142
        err_mat=[]
# ddd/Package_P_Matrix_Comparisons.m:146
    
    Grad_mat[abs(Grad_mat) > 10 ** 5]=0
# ddd/Package_P_Matrix_Comparisons.m:155
    #disp(error)
    
    ## time keeping
    stp_time=copy(toc)
# ddd/Package_P_Matrix_Comparisons.m:160
    if stp_time > 5:
        disp(strcat('slow eval:',num2str(stp_time)))
    
    return error,Grad_mat,err_mat
    
if __name__ == '__main__':
    pass
    