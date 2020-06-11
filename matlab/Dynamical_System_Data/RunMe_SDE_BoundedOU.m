clear all
close all

dt = 2^-9;

l_lim = 0;
u_lim = 5;

start_time = 0;
end_time = 100;
n_gridpoints = ((end_time - start_time)/dt) + 1; 
time_vec = linspace(start_time,end_time,n_gridpoints);

syms x y
f = symfun( 0.5*exp(-2*((x-0.5)^2 + (y-0.5)^2)) + 1*exp(-1*((x-4)^2 + (y-2)^2)) + 1*exp(-1*((x-2)^2 + (y-4)^2)) -1*( ( 0.5*(x^2 + y^2) - x*y ) - 0.1*(x^2+y^2) )^2  , [x y]);

dfdx = diff(f,x);
dfdy = diff(f,y);
F_pot = matlabFunction(f);
V_x = matlabFunction(dfdx);
V_y = matlabFunction(dfdy);

ln = linspace(l_lim,u_lim,101);
[x_grid,y_grid] = ndgrid(ln,ln);
%%
c_lim=1.5;
figure(1), hold on
[a,b] = contourf(x_grid,y_grid,-F_pot(x_grid,y_grid), linspace(-1,c_lim,21)); alpha(0.25)
set(b,'LineColor','none')
caxis([-1,c_lim])
%set(gca,'EdgeAlpha',0,'FaceAlpha',0.7);
%quiver(x_grid,y_grid,V_x(x_grid,y_grid),V_y(x_grid,y_grid),'k')
   set(gca,'fontsize',17)
axis tight square
axis([l_lim,u_lim,l_lim,u_lim])
tightfig();
h = colorbar;
h.Location = 'north';
h.Color = [0,0,0];
h.FontSize = 15;
colormap hot
box on
%%
Diff_const = 0.25;% Diff_const = 0.075;
sigma_const_1 = sqrt(2*Diff_const);
sigma_const_2 = sqrt(2*Diff_const);
%
Num_pts = 1000;
sol_X_1 = 0.5*randn(Num_pts,n_gridpoints) + 0.5;
sol_X_2 = 0.5*randn(Num_pts,n_gridpoints) + 0.5;
ind_1 = sol_X_1(:,1) < 0;
sol_X_1(ind_1,1) = - sol_X_1(ind_1,1);
ind_2 = sol_X_2(:,1) < 0;
sol_X_2(ind_2,1) = - sol_X_2(ind_2,1);

%
for ii = 1:(n_gridpoints-1)
    rand_noise_1 = sigma_const_1.*sqrt(dt).*randn(Num_pts,1);
    rand_noise_2 = sigma_const_2.*sqrt(dt).*randn(Num_pts,1);
    sol_X_1(:,ii+1) = sol_X_1(:,ii) + V_x(sol_X_1(:,ii),sol_X_2(:,ii)).*dt + rand_noise_1;
    sol_X_2(:,ii+1) = sol_X_2(:,ii) + V_y(sol_X_1(:,ii),sol_X_2(:,ii)).*dt + rand_noise_2;
    
     ind_1 = sol_X_1(:,ii+1) < 0;
     sol_X_1(ind_1,ii+1) = - sol_X_1(ind_1,ii+1);
%         
     ind_2 = sol_X_2(:,ii+1) < 0;
     sol_X_2(ind_2,ii+1) = - sol_X_2(ind_2,ii+1);
    
end
%


%%
ind = [1,2,3];
figure()
plot(sol_X_1(ind,1:find(time_vec == 89))',sol_X_2(ind,1:find(time_vec == 89))');
   set(gca,'fontsize',17)
axis tight square

axis([l_lim,u_lim,l_lim,u_lim])
   set(gca,'XTick',l_lim:u_lim)
   set(gca,'YTick',l_lim:u_lim)
   tightfig();
%%
figure(), hold on
plot(std(sol_X_1),std(sol_X_2),'r')
plot(mean(sol_X_1),mean(sol_X_2),'k')
plot(linspace(0,4,51),linspace(0,4,51))

%%


scale_fac = 2^7; %Option of down-sampling time.
save_traj_1 = sol_X_1(:,1:scale_fac:n_gridpoints);
save_traj_2 = sol_X_2(:,1:scale_fac:n_gridpoints);
save_dt = dt*scale_fac;
%%

%save('SDE_2d_Simulated_Switch', 'save_traj_1','save_traj_2', 'save_dt' )