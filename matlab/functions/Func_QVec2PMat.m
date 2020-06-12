function P_Mat = Func_QVec2PMat(Q_Vec, Model)

Start_Vec = Q_Vec(end-Model.num_basis+1:end);
Vec_form = Q_Vec(1:end-Model.num_basis);
Q_Mat = double( Model.Indc );
Q_Mat( Model.Indc ) = Vec_form;
P_Mat = Model.MassMat\Q_Mat;
P_Mat(:,end+1) = Start_Vec;
P_Mat = full(P_Mat);

end
