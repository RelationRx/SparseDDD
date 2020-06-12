%% Preprocess SDE and put into file format for MATLAB version of Sparse DDD + Create CSV files for python version.
clear all
close all
clc
cd ../
cd simulate_dynamical_system/
load('SDE_2d_Simulated_Switch.mat')
cd ../
cd preprocess/
%% Define observation points.
t_map = 0:save_dt:100;
Global.time_pts = [0,1,2,3,5,8,13,21,34,55,89];
Global.num_tpts = numel(Global.time_pts);
Global.Num_pts_each_time = 1000*ones(size(Global.time_pts));
[~,idxs_of_interest] = intersect(t_map,Global.time_pts);
Global.Data_Dim = 2;
%% Reduce time series to observation points only
X1_sol_t_pts = save_traj_1(:,idxs_of_interest);
X2_sol_t_pts = save_traj_2(:,idxs_of_interest);
%% Create Global variable to save both individual time points and concatenated dataset for ease.
Global.sample_pts = cell(numel(Global.time_pts),1);
Global.all_pts_MixModel = zeros(sum(Global.Num_pts_each_time),Global.Data_Dim);
Global.running_indx_sum = [0,cumsum(Global.Num_pts_each_time)];
for uu = 1:numel(Global.time_pts)
    Global.sample_pts{uu} = [X1_sol_t_pts(:,uu),X2_sol_t_pts(:,uu)];
    Global.all_pts_MixModel(Global.running_indx_sum(uu)+1:Global.running_indx_sum(uu+1),:) = Global.sample_pts{uu};
    disp( sum(Global.sample_pts{uu}(:,1) > Global.sample_pts{uu}(:,2)) )
end
clear uu
%% Save Global variable file.
save('SDE_ProcessedData.mat','Global')
%% Save CSVs for python code.
cd ../
cd raw_csv_files/
for ii = 1:11
    writematrix(Global.sample_pts{ii},strcat('SDE_observed_at_time_',num2str(Global.time_pts(ii)),'.csv') );
end