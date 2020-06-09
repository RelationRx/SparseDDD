# Generated with SMOP  0.41-beta
from libsmop import *
# ddd/Func_CalculateGMatrix.m

    
@function
def Func_CalculateGMatrix(PF_mat=None,tmp_middle_mat=None,t_gap=None,normalisation_fac=None,*args,**kwargs):
    varargin = Func_CalculateGMatrix.varargin
    nargin = Func_CalculateGMatrix.nargin

    PFt=PF_mat.T
# ddd/Func_CalculateGMatrix.m:2
    tol=1e-10
# ddd/Func_CalculateGMatrix.m:3
    S_k0=copy(tmp_middle_mat)
# ddd/Func_CalculateGMatrix.m:4
    S_k1=dot(PFt,S_k0) + dot(S_k0,PFt)
# ddd/Func_CalculateGMatrix.m:5
    # Zeroth team
    scale_fac_tmp=t_gap / 1
# ddd/Func_CalculateGMatrix.m:7
    tmp_Grad_mat=dot(S_k0,scale_fac_tmp)
# ddd/Func_CalculateGMatrix.m:8
    # First term
    scale_fac_tmp=dot(t_gap,scale_fac_tmp) / 2
# ddd/Func_CalculateGMatrix.m:10
    tmp_Grad_mat=tmp_Grad_mat + dot(S_k1,scale_fac_tmp)
# ddd/Func_CalculateGMatrix.m:11
    trigger=0
# ddd/Func_CalculateGMatrix.m:12
    aa=2
# ddd/Func_CalculateGMatrix.m:13
    # Second term onwards.
    while trigger == 0:

        S_k2=dot(PFt,S_k1) + dot(S_k1,PFt) - dot(dot(PFt,S_k0),PFt)
# ddd/Func_CalculateGMatrix.m:16
        scale_fac_tmp=dot(t_gap,scale_fac_tmp) / (aa + 1)
# ddd/Func_CalculateGMatrix.m:17
        Mat_to_add=dot(S_k2,scale_fac_tmp)
# ddd/Func_CalculateGMatrix.m:18
        if any(logical_not(isfinite(ravel(Mat_to_add)))):
            break
        tmp_Grad_mat=tmp_Grad_mat + Mat_to_add
# ddd/Func_CalculateGMatrix.m:24
        S_k0=copy(S_k1)
# ddd/Func_CalculateGMatrix.m:26
        S_k1=copy(S_k2)
# ddd/Func_CalculateGMatrix.m:27
        if max(ravel(Mat_to_add)) < dot(tol,normalisation_fac):
            trigger=copy(aa)
# ddd/Func_CalculateGMatrix.m:29
        else:
            aa=aa + 1
# ddd/Func_CalculateGMatrix.m:31

    
    G=dot(2,tmp_Grad_mat) / normalisation_fac
# ddd/Func_CalculateGMatrix.m:34
    return G
    
if __name__ == '__main__':
    pass
    