import math
import numpy


def multiply(x, y):
    return x * y


def mc_integration(ii=None, jj=None, Model=None, Global=None, X0=None, *args, **kwargs):
    max_rad = max(Model.Basis_eps(ii), Model.Basis_eps(jj))
    min_rad = min(Model.Basis_eps(ii), Model.Basis_eps(jj))
    offset = Model.Basis_Loc[jj] - Model.Basis_Loc[ii]
    radnorm = numpy.linalg.norm(offset, 2)
    if radnorm < (max_rad + min_rad):
        x_samp = ((2 * max_rad) * X0) - max_rad
        if Model.Basis_eps(ii) > Model.Basis_eps(jj):
            tmp = multiply(Model.basis_nd(x_samp, Global.Data_Dim, Model.Basis_eps(ii)),
                                        Model.basis_nd(x_samp - offset, Global.Data_Dim, Model.Basis_eps(jj)))
        else:
            tmp = multiply(Model.basis_nd(x_samp + offset, Global.Data_Dim, Model.Basis_eps(ii)),
                                        Model.basis_nd(x_samp, Global.Data_Dim, Model.Basis_eps(jj)))
        value = (((2 * max_rad) ** Global.Data_Dim) * sum(tmp)) / len(tmp)
    else:
        value = 0
    return value
