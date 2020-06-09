# Generated with SMOP  0.41-beta
from libsmop import *
# ddd/Package_Func_SingleTimePointError.m

    
@function
def Package_Func_SingleTimePointError(PF_Mat=None,Coeff_vec=None,Start_Vec=None,t_gap=None,Model=None,*args,**kwargs):
    varargin = Package_Func_SingleTimePointError.varargin
    nargin = Package_Func_SingleTimePointError.nargin

    if Model.norm_fac == 'Relative':
        normalisation_fac=(dot(dot(Coeff_vec.T,Model.MassMat),Coeff_vec))
# ddd/Package_Func_SingleTimePointError.m:3
    else:
        normalisation_fac=1
# ddd/Package_Func_SingleTimePointError.m:5
    
    #disp(normalisation_fac)
    mat_G=expm(dot(t_gap,PF_Mat))
# ddd/Package_Func_SingleTimePointError.m:8
    result_vec=dot(mat_G,Start_Vec) - Coeff_vec
# ddd/Package_Func_SingleTimePointError.m:9
    error=(dot(dot(result_vec.T,Model.MassMat),result_vec)) / normalisation_fac
# ddd/Package_Func_SingleTimePointError.m:10
    if isnan(error):
        error=copy(nan)
# ddd/Package_Func_SingleTimePointError.m:12
    
    if isnan(error):
        Grad_mat=nan(size(PF_Mat,1),size(PF_Mat,1) + 1)
# ddd/Package_Func_SingleTimePointError.m:15
    else:
        tmp_middle_mat=dot(dot(Model.MassMat,result_vec),Start_Vec.T)
# ddd/Package_Func_SingleTimePointError.m:17
        G_tmp=Func_CalculateGMatrix(PF_Mat,tmp_middle_mat,t_gap,normalisation_fac)
# ddd/Package_Func_SingleTimePointError.m:18
        G_tmp_2=(dot(dot(dot(2,(mat_G.T)),Model.MassMat),result_vec)) / normalisation_fac
# ddd/Package_Func_SingleTimePointError.m:19
        Grad_mat=concat([G_tmp,G_tmp_2])
# ddd/Package_Func_SingleTimePointError.m:20
    
    return error,Grad_mat
    
if __name__ == '__main__':
    pass
    