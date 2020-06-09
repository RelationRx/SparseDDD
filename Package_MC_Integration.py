# Generated with SMOP  0.41-beta
from libsmop import *
# ddd/Package_MC_Integration.m

    
@function
def Package_MC_Integration(ii=None,jj=None,Model=None,Global=None,X0=None,*args,**kwargs):
    varargin = Package_MC_Integration.varargin
    nargin = Package_MC_Integration.nargin

    max_rad=max(Model.Basis_eps(ii),Model.Basis_eps(jj))
# ddd/Package_MC_Integration.m:3
    min_rad=min(Model.Basis_eps(ii),Model.Basis_eps(jj))
# ddd/Package_MC_Integration.m:4
    offset=Model.Basis_Loc(jj,arange()) - Model.Basis_Loc(ii,arange())
# ddd/Package_MC_Integration.m:5
    radnorm=norm(offset,2)
# ddd/Package_MC_Integration.m:6
    if radnorm < (max_rad + min_rad):
        if Model.Basis_eps(ii) > Model.Basis_eps(jj):
            f=lambda x=None: multiply(Model.basis_nd(x,Global.Data_Dim,Model.Basis_eps(ii)),Model.basis_nd(x - offset,Global.Data_Dim,Model.Basis_eps(jj)))
# ddd/Package_MC_Integration.m:10
        else:
            f=lambda x=None: multiply(Model.basis_nd(x + offset,Global.Data_Dim,Model.Basis_eps(ii)),Model.basis_nd(x,Global.Data_Dim,Model.Basis_eps(jj)))
# ddd/Package_MC_Integration.m:12
        X_samp=(dot(dot(2,max_rad),X0)) - max_rad
# ddd/Package_MC_Integration.m:15
        tmp=f(X_samp)
# ddd/Package_MC_Integration.m:16
        tmp=dot(((dot(2,max_rad)) ** Global.Data_Dim),sum(tmp)) / length(tmp)
# ddd/Package_MC_Integration.m:17
        Value=copy(tmp)
# ddd/Package_MC_Integration.m:18
    else:
        Value=0
# ddd/Package_MC_Integration.m:20
    
    return Value
    
if __name__ == '__main__':
    pass
    