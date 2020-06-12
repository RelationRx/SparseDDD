function [error,Grad_mat] = Package_Func_SingleTimePointError(PF_Mat,Coeff_vec,Start_Vec, t_gap, Model)
if Model.norm_fac == "Relative"
    normalisation_fac = (Coeff_vec'*Model.MassMat*Coeff_vec);
else
    normalisation_fac = 1;
end
%disp(normalisation_fac)
mat_G = expm(t_gap*PF_Mat);
result_vec = mat_G*Start_Vec - Coeff_vec;
error = (result_vec'*Model.MassMat*result_vec)/normalisation_fac;
if isnan(error)
    error = nan;
end
if isnan(error)
    Grad_mat = nan(size(PF_Mat,1),size(PF_Mat,1) + 1);
else
    tmp_middle_mat = Model.MassMat*result_vec*Start_Vec';
    G_tmp = Func_CalculateGMatrix(PF_Mat,tmp_middle_mat,t_gap,normalisation_fac);
    G_tmp_2 = (2*(mat_G')*Model.MassMat*result_vec)/normalisation_fac;
    Grad_mat = [G_tmp,G_tmp_2];
end
end
