
dt=2 ** - 9
l_lim=0
u_lim=5
start_time=0
end_time=100
n_gridpoints=((end_time - start_time) / dt) + 1
time_vec=linspace(start_time,end_time,n_gridpoints)
syms('x','y')
f=symfun(dot(0.5,exp(dot(- 2,((x - 0.5) ** 2 + (y - 0.5) ** 2)))) + dot(1,exp(dot(- 1,((x - 4) ** 2 + (y - 2) ** 2)))) + dot(1,exp(dot(- 1,((x - 2) ** 2 + (y - 4) ** 2)))) - dot(1,((dot(0.5,(x ** 2 + y ** 2)) - dot(x,y)) - dot(0.1,(x ** 2 + y ** 2))) ** 2),concat([x,y]))
dfdx=diff(f,x)
dfdy=diff(f,y)
F_pot=matlabFunction(f)
V_x=matlabFunction(dfdx)
V_y=matlabFunction(dfdy)
ln=linspace(l_lim,u_lim,101)
x_grid,y_grid=ndgrid(ln,ln,nargout=2)
##
c_lim=1.5
a,b=contourf(x_grid,y_grid,- F_pot(x_grid,y_grid),linspace(- 1,c_lim,21),nargout=2)
alpha(0.25)
h=copy(colorbar)
h.Location = copy('north')
h.Color = copy(concat([0,0,0]))
h.FontSize = copy(15)
colormap('hot')
box('on')
##
Diff_const=0.25
sigma_const_1=sqrt(dot(2,Diff_const))
sigma_const_2=sqrt(dot(2,Diff_const))

Num_pts=1000
sol_X_1=dot(0.5,randn(Num_pts,n_gridpoints)) + 0.5
sol_X_2=dot(0.5,randn(Num_pts,n_gridpoints)) + 0.5
ind_1=sol_X_1(arange(),1) < 0
sol_X_1[ind_1,1]=- sol_X_1(ind_1,1)
ind_2=sol_X_2(arange(),1) < 0
sol_X_2[ind_2,1]=- sol_X_2(ind_2,1)

for ii in arange(1,(n_gridpoints - 1)).reshape(-1):
    rand_noise_1=multiply(multiply(sigma_const_1,sqrt(dt)),randn(Num_pts,1))
    rand_noise_2=multiply(multiply(sigma_const_2,sqrt(dt)),randn(Num_pts,1))
    sol_X_1[arange(),ii + 1]=sol_X_1(arange(),ii) + multiply(V_x(sol_X_1(arange(),ii),sol_X_2(arange(),ii)),dt) + rand_noise_1
    sol_X_2[arange(),ii + 1]=sol_X_2(arange(),ii) + multiply(V_y(sol_X_1(arange(),ii),sol_X_2(arange(),ii)),dt) + rand_noise_2
    ind_1=sol_X_1(arange(),ii + 1) < 0
    sol_X_1[ind_1,ii + 1]=- sol_X_1(ind_1,ii + 1)
    #
    ind_2=sol_X_2(arange(),ii + 1) < 0
    sol_X_2[ind_2,ii + 1]=- sol_X_2(ind_2,ii + 1)


##
ind=concat([1,2,3])
figure()
plot(sol_X_1(ind,arange(1,find(time_vec == 89))).T,sol_X_2(ind,arange(1,find(time_vec == 89))).T)
set(gca,'fontsize',17)
axis('tight','square')
axis(concat([l_lim,u_lim,l_lim,u_lim]))
set(gca,'XTick',arange(l_lim,u_lim))
set(gca,'YTick',arange(l_lim,u_lim))
tightfig()
##
figure()
hold('on')
plot(std(sol_X_1),std(sol_X_2),'r')
plot(mean(sol_X_1),mean(sol_X_2),'k')
plot(linspace(0,4,51),linspace(0,4,51))
##

scale_fac=2 ** 7

save_traj_1=sol_X_1(arange(),arange(1,n_gridpoints,scale_fac))
save_traj_2=sol_X_2(arange(),arange(1,n_gridpoints,scale_fac))
save_dt=dot(dt,scale_fac)
##

#save('SDE_2d_Simulated_Switch', 'save_traj_1','save_traj_2', 'save_dt' )