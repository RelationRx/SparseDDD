clear all
close all
addpath('Dynamical_System_Data/')
load('SDE_2d_Simulated_Switch.mat')

%%
t_map = 0:save_dt:100;
Global.time_pts = [0,1,2,3,5,8,13,21,34,55,89];
Global.num_tpts = numel(Global.time_pts);
Global.Num_pts_each_time = 1000*ones(size(Global.time_pts));
[~,idxs_of_interest] = intersect(t_map,Global.time_pts);
%%
X1_sol_t_pts = save_traj_1(:,idxs_of_interest);
X2_sol_t_pts = save_traj_2(:,idxs_of_interest);
%%
Global.Data_Dim = 2;
%%
Global.sample_pts = cell(numel(Global.time_pts),1);
Global.all_pts_MixModel = zeros(sum(Global.Num_pts_each_time),Global.Data_Dim);
Global.running_indx_sum = [0,cumsum(Global.Num_pts_each_time)];
for uu = 1:numel(Global.time_pts)
    Global.sample_pts{uu} = [X1_sol_t_pts(:,uu),X2_sol_t_pts(:,uu)];
    Global.all_pts_MixModel(Global.running_indx_sum(uu)+1:Global.running_indx_sum(uu+1),:) = Global.sample_pts{uu};
    disp( sum(Global.sample_pts{uu}(:,1) > Global.sample_pts{uu}(:,2)) )
end
clear uu

%%
save('SDE_ProcessedData.mat','Global')
