# Generated with SMOP  0.41-beta
from libsmop import *
# ddd/Package_Q_Matrix_Comparisons.m

    
@function
def Package_Q_Matrix_Comparisons(Q_Mat=None,Model=None,Global=None,*args,**kwargs):
    varargin = Package_Q_Matrix_Comparisons.varargin
    nargin = Package_Q_Matrix_Comparisons.nargin

    tic
    Start_Vec=Q_Mat(arange(end() - Model.num_basis + 1,end()))
# ddd/Package_Q_Matrix_Comparisons.m:4
    Vec_form=Q_Mat(arange(1,end() - Model.num_basis))
# ddd/Package_Q_Matrix_Comparisons.m:5
    Q_Mat=double(Model.Indc)
# ddd/Package_Q_Matrix_Comparisons.m:6
    Q_Mat[Model.Indc]=Vec_form
# ddd/Package_Q_Matrix_Comparisons.m:7
    PF_Mat=numpy.linalg.solve(Model.MassMat,Q_Mat)
# ddd/Package_Q_Matrix_Comparisons.m:8
    ## Full time
    if Model.ShowComparison == 'FullSnapshot':
        err_mat=zeros(numel(Global.time_pts),1)
# ddd/Package_Q_Matrix_Comparisons.m:14
        Grad_mat=zeros(size(Q_Mat,1),size(Q_Mat,1) + 1)
# ddd/Package_Q_Matrix_Comparisons.m:15
        for vv in arange(1,numel(Global.time_pts)).reshape(-1):
            t_gap=(Global.time_pts(vv) - Global.time_pts(1))
# ddd/Package_Q_Matrix_Comparisons.m:17
            Coeff_vec=Model.Coeff_vecs[vv]
# ddd/Package_Q_Matrix_Comparisons.m:18
            error_tmp,Grad_mat_tmp=Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Start_Vec,t_gap,Model,nargout=2)
# ddd/Package_Q_Matrix_Comparisons.m:19
            err_mat[vv]=error_tmp
# ddd/Package_Q_Matrix_Comparisons.m:20
            Grad_mat=Grad_mat + Grad_mat_tmp
# ddd/Package_Q_Matrix_Comparisons.m:21
        error=sum(err_mat) / Global.num_tpts
# ddd/Package_Q_Matrix_Comparisons.m:24
        Grad_mat=Grad_mat / Global.num_tpts
# ddd/Package_Q_Matrix_Comparisons.m:25
        Grad_mat[arange(),arange(1,Model.num_basis)]=numpy.linalg.solve(Model.MassMat,Grad_mat(arange(),arange(1,Model.num_basis)))
# ddd/Package_Q_Matrix_Comparisons.m:27
        Grad_mat=Grad_mat(logical(concat([Model.Indc,ones(Model.num_basis,1)])))
# ddd/Package_Q_Matrix_Comparisons.m:28
    
    ## Evaluate first and last point only!
    if Model.ShowComparison == 'FirstLastPoint':
        err_mat=zeros(2,1)
# ddd/Package_Q_Matrix_Comparisons.m:36
        Grad_mat=zeros(size(Q_Mat,1),size(Q_Mat,1) + 1)
# ddd/Package_Q_Matrix_Comparisons.m:37
        for vv in Model.pts_of_int.reshape(-1):
            t_gap=(Global.time_pts(vv) - Global.time_pts(1))
# ddd/Package_Q_Matrix_Comparisons.m:39
            Coeff_vec=Model.Coeff_vecs[vv]
# ddd/Package_Q_Matrix_Comparisons.m:40
            error_tmp,Grad_mat_tmp=Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Start_Vec,t_gap,Model,nargout=2)
# ddd/Package_Q_Matrix_Comparisons.m:41
            err_mat[vv]=error_tmp
# ddd/Package_Q_Matrix_Comparisons.m:42
            Grad_mat=Grad_mat + Grad_mat_tmp
# ddd/Package_Q_Matrix_Comparisons.m:43
        err_mat=err_mat(Model.pts_of_int)
# ddd/Package_Q_Matrix_Comparisons.m:45
        error=sum(err_mat) / 2
# ddd/Package_Q_Matrix_Comparisons.m:46
        Grad_mat=Grad_mat / 2
# ddd/Package_Q_Matrix_Comparisons.m:47
        Grad_mat[arange(),arange(1,Model.num_basis)]=numpy.linalg.solve(Model.MassMat,Grad_mat(arange(),arange(1,Model.num_basis)))
# ddd/Package_Q_Matrix_Comparisons.m:48
        Grad_mat=Grad_mat(logical(concat([Model.Indc,ones(Model.num_basis,1)])))
# ddd/Package_Q_Matrix_Comparisons.m:49
    
    ## Time Series
    if Model.ShowComparison == 'FullTimeSeries':
        err_mat=zeros(concat([size(Model.TimeSeries,1) - 1,Model.n_ts_samples]))
# ddd/Package_Q_Matrix_Comparisons.m:59
        Grad_mat=zeros(size(Q_Mat,1),size(Q_Mat,1) + 1)
# ddd/Package_Q_Matrix_Comparisons.m:60
        tot_valid=0
# ddd/Package_Q_Matrix_Comparisons.m:61
        for ww in arange(1,Model.n_ts_samples).reshape(-1):
            #for ww = randsample(1000,Model.n_ts_samples)'
            err_vec=zeros(1,numel(Global.time_pts) - 1)
# ddd/Package_Q_Matrix_Comparisons.m:64
            Grad_mat_cache=zeros(size(Q_Mat,1),size(Q_Mat,1) + 1)
# ddd/Package_Q_Matrix_Comparisons.m:65
            for vv in arange(2,numel(Global.time_pts)).reshape(-1):
                t_gap=(Global.time_pts(vv) - Global.time_pts(1))
# ddd/Package_Q_Matrix_Comparisons.m:67
                Coeff_vec=Model.TimeSeries[vv,ww]
# ddd/Package_Q_Matrix_Comparisons.m:68
                error_tmp,Grad_mat_tmp=Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Model.TimeSeries[1,ww],t_gap,Model,nargout=2)
# ddd/Package_Q_Matrix_Comparisons.m:69
                err_vec[vv - 1]=error_tmp
# ddd/Package_Q_Matrix_Comparisons.m:70
                Grad_mat_cache=Grad_mat_cache + Grad_mat_tmp
# ddd/Package_Q_Matrix_Comparisons.m:71
            if sum(isnan(err_vec)) == 0:
                err_mat[arange(),ww]=err_vec
# ddd/Package_Q_Matrix_Comparisons.m:74
                Grad_mat=Grad_mat + Grad_mat_cache
# ddd/Package_Q_Matrix_Comparisons.m:75
                tot_valid=tot_valid + 1
# ddd/Package_Q_Matrix_Comparisons.m:76
        err_mat[arange(),sum(err_mat) == 0]=[]
# ddd/Package_Q_Matrix_Comparisons.m:80
        error=sum(sum(err_mat)) / (dot((Global.num_tpts - 1),tot_valid))
# ddd/Package_Q_Matrix_Comparisons.m:82
        Grad_mat=Grad_mat / (dot((Global.num_tpts - 1),tot_valid))
# ddd/Package_Q_Matrix_Comparisons.m:83
        Grad_mat[arange(),arange(1,Model.num_basis)]=numpy.linalg.solve(Model.MassMat,Grad_mat(arange(),arange(1,Model.num_basis)))
# ddd/Package_Q_Matrix_Comparisons.m:87
        Grad_mat[arange(),end()]=0
# ddd/Package_Q_Matrix_Comparisons.m:88
        Grad_mat=Grad_mat(logical(concat([Model.Indc,ones(Model.num_basis,1)])))
# ddd/Package_Q_Matrix_Comparisons.m:90
    
    ## Comparison
    if Model.ShowComparison == 'Composite':
        if Model.Lambda != 0:
            err_mat_TS=zeros(concat([size(Model.TimeSeries,1) - 1,Model.n_ts_samples]))
# ddd/Package_Q_Matrix_Comparisons.m:99
            Grad_mat_TS=zeros(size(Q_Mat,1),size(Q_Mat,1) + 1)
# ddd/Package_Q_Matrix_Comparisons.m:100
            tot_valid=0
# ddd/Package_Q_Matrix_Comparisons.m:101
            for ww in arange(1,Model.n_ts_samples).reshape(-1):
                err_vec=zeros(1,numel(Global.time_pts) - 1)
# ddd/Package_Q_Matrix_Comparisons.m:103
                Grad_mat_cache=zeros(size(Q_Mat,1),size(Q_Mat,1) + 1)
# ddd/Package_Q_Matrix_Comparisons.m:104
                for vv in arange(2,numel(Global.time_pts)).reshape(-1):
                    t_gap=(Global.time_pts(vv) - Global.time_pts(1))
# ddd/Package_Q_Matrix_Comparisons.m:106
                    Coeff_vec=Model.TimeSeries[vv,ww]
# ddd/Package_Q_Matrix_Comparisons.m:107
                    error_tmp,Grad_mat_tmp=Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Model.TimeSeries[1,ww],t_gap,Model,nargout=2)
# ddd/Package_Q_Matrix_Comparisons.m:108
                    err_vec[vv - 1]=error_tmp
# ddd/Package_Q_Matrix_Comparisons.m:109
                    Grad_mat_cache=Grad_mat_cache + Grad_mat_tmp
# ddd/Package_Q_Matrix_Comparisons.m:110
                if sum(isnan(err_vec)) == 0:
                    err_mat_TS[arange(),ww]=err_vec
# ddd/Package_Q_Matrix_Comparisons.m:113
                    Grad_mat_TS=Grad_mat_TS + Grad_mat_cache
# ddd/Package_Q_Matrix_Comparisons.m:114
                    tot_valid=tot_valid + 1
# ddd/Package_Q_Matrix_Comparisons.m:115
            err_mat_TS[arange(),sum(err_mat_TS) == 0]=[]
# ddd/Package_Q_Matrix_Comparisons.m:118
            error_TS=sum(sum(err_mat_TS)) / (dot((Global.num_tpts - 1),tot_valid))
# ddd/Package_Q_Matrix_Comparisons.m:120
            Grad_mat_TS=Grad_mat_TS / (dot((Global.num_tpts - 1),tot_valid))
# ddd/Package_Q_Matrix_Comparisons.m:121
            Grad_mat_TS[arange(),arange(1,Model.num_basis)]=numpy.linalg.solve(Model.MassMat,Grad_mat_TS(arange(),arange(1,Model.num_basis)))
# ddd/Package_Q_Matrix_Comparisons.m:123
            Grad_mat_TS[arange(),end()]=0
# ddd/Package_Q_Matrix_Comparisons.m:124
        # First last Snapshot
        err_mat_FL=zeros(2,1)
# ddd/Package_Q_Matrix_Comparisons.m:130
        Grad_mat_FL=zeros(size(Q_Mat,1),size(Q_Mat,1) + 1)
# ddd/Package_Q_Matrix_Comparisons.m:131
        for vv in Model.pts_of_int.reshape(-1):
            t_gap=(Global.time_pts(vv) - Global.time_pts(1))
# ddd/Package_Q_Matrix_Comparisons.m:133
            Coeff_vec=Model.Coeff_vecs[vv]
# ddd/Package_Q_Matrix_Comparisons.m:134
            error_tmp,Grad_mat_tmp=Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Start_Vec,t_gap,Model,nargout=2)
# ddd/Package_Q_Matrix_Comparisons.m:135
            err_mat_FL[vv]=error_tmp
# ddd/Package_Q_Matrix_Comparisons.m:136
            Grad_mat_FL=Grad_mat_FL + Grad_mat_tmp
# ddd/Package_Q_Matrix_Comparisons.m:137
        err_mat_FL=err_mat_FL(Model.pts_of_int)
# ddd/Package_Q_Matrix_Comparisons.m:139
        error_FL=sum(err_mat_FL) / 2
# ddd/Package_Q_Matrix_Comparisons.m:140
        Grad_mat_FL=Grad_mat_FL / 2
# ddd/Package_Q_Matrix_Comparisons.m:141
        Grad_mat_FL[arange(),arange(1,Model.num_basis)]=numpy.linalg.solve(Model.MassMat,Grad_mat_FL(arange(),arange(1,Model.num_basis)))
# ddd/Package_Q_Matrix_Comparisons.m:142
        if Model.Lambda != 0:
            error=dot(Model.Lambda,error_TS) + dot((1 - Model.Lambda),error_FL)
# ddd/Package_Q_Matrix_Comparisons.m:147
            Grad_mat=dot(Model.Lambda,Grad_mat_TS) + dot((1 - Model.Lambda),Grad_mat_FL)
# ddd/Package_Q_Matrix_Comparisons.m:148
            disp(concat([error_TS,error_FL,error]))
        else:
            error=copy(error_FL)
# ddd/Package_Q_Matrix_Comparisons.m:151
            Grad_mat=copy(Grad_mat_FL)
# ddd/Package_Q_Matrix_Comparisons.m:152
        err_mat=[]
# ddd/Package_Q_Matrix_Comparisons.m:156
        Grad_mat=Grad_mat(logical(concat([Model.Indc,ones(Model.num_basis,1)])))
# ddd/Package_Q_Matrix_Comparisons.m:161
    
    Grad_mat[abs(Grad_mat) > 10 ** 5]=0
# ddd/Package_Q_Matrix_Comparisons.m:166
    ## time keeping
    stp_time=copy(toc)
# ddd/Package_Q_Matrix_Comparisons.m:171
    if stp_time > 5:
        disp(strcat('slow eval:',num2str(stp_time)))
    
    return error,Grad_mat,err_mat
    
if __name__ == '__main__':
    pass
    