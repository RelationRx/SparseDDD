# Generated with SMOP  0.41-beta
from libsmop import *
# ddd/Package_Func_QVec2PMat.m

    
@function
def Package_Func_QVec2PMat(Q_Vec=None,Model=None,*args,**kwargs):
    varargin = Package_Func_QVec2PMat.varargin
    nargin = Package_Func_QVec2PMat.nargin

    Start_Vec=Q_Vec(arange(end() - Model.num_basis + 1,end()))
# ddd/Package_Func_QVec2PMat.m:3
    Vec_form=Q_Vec(arange(1,end() - Model.num_basis))
# ddd/Package_Func_QVec2PMat.m:4
    Q_Mat=double(Model.Indc)
# ddd/Package_Func_QVec2PMat.m:5
    Q_Mat[Model.Indc]=Vec_form
# ddd/Package_Func_QVec2PMat.m:6
    P_Mat=numpy.linalg.solve(Model.MassMat,Q_Mat)
# ddd/Package_Func_QVec2PMat.m:7
    P_Mat[arange(),end() + 1]=Start_Vec
# ddd/Package_Func_QVec2PMat.m:8
    P_Mat=full(P_Mat)
# ddd/Package_Func_QVec2PMat.m:9
    return P_Mat
    
if __name__ == '__main__':
    pass
    