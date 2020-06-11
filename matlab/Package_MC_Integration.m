function Value = Package_MC_Integration(ii,jj,Model, Global, X0)

max_rad = max( Model.Basis_eps(ii), Model.Basis_eps(jj) );
min_rad = min( Model.Basis_eps(ii), Model.Basis_eps(jj) );
offset = Model.Basis_Loc(jj,:) - Model.Basis_Loc(ii,:);
radnorm = norm( offset, 2);
if radnorm < ( max_rad + min_rad )
    
if Model.Basis_eps(ii) > Model.Basis_eps(jj)     
    f = @(x) Model.basis_nd(x ,Global.Data_Dim, Model.Basis_eps(ii) ).*Model.basis_nd(x - offset,Global.Data_Dim, Model.Basis_eps(jj) );
else
    f = @(x) Model.basis_nd(x + offset,Global.Data_Dim, Model.Basis_eps(ii) ).*Model.basis_nd(x,Global.Data_Dim, Model.Basis_eps(jj) );
end
    
X_samp = (2*max_rad*X0) - max_rad;
tmp = f(X_samp);
tmp =  ((2*max_rad)^Global.Data_Dim) * sum(tmp) / length(tmp);
Value = tmp;
else
    Value = 0;
end
end

