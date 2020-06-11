function G = Func_CalculateGMatrix(PF_mat,tmp_middle_mat,t_gap,normalisation_fac)
PFt = PF_mat';
tol = 1e-10;
S_k0 = tmp_middle_mat;
S_k1 = PFt*S_k0 + S_k0*PFt;
% Zeroth team
scale_fac_tmp = t_gap/1;
tmp_Grad_mat = S_k0*scale_fac_tmp; 
% First term
scale_fac_tmp = t_gap*scale_fac_tmp/2;
tmp_Grad_mat = tmp_Grad_mat + S_k1*scale_fac_tmp;
trigger = 0;
aa = 2;
% Second term onwards.
while trigger == 0 
    S_k2 = PFt*S_k1 + S_k1*PFt - PFt*S_k0*PFt;
    scale_fac_tmp = t_gap*scale_fac_tmp/(aa+1);
    Mat_to_add = S_k2*scale_fac_tmp;
    
    if any(~isfinite(Mat_to_add(:)))
        break
    end
    
    tmp_Grad_mat = tmp_Grad_mat + Mat_to_add;
    
    S_k0 = S_k1;
    S_k1 = S_k2;
    if max(Mat_to_add(:)) < tol*normalisation_fac
        trigger = aa;
    else
        aa = aa + 1;
    end    
end
G = 2*tmp_Grad_mat/normalisation_fac;
end