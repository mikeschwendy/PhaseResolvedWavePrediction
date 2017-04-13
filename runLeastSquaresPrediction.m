function [z_target_pred,z_target_truth,t_pred] = runLeastSquaresPrediction(...
    x_target,y_target,dt,x_array,y_array,z_array,...
    k,theta,reg_factor,T_meas,T_pred,T_delay,overlap,Y)


% reshape/initialize arrays
nx_target = length(x_target);
ny_target = length(y_target);
[x_target,y_target] = meshgrid(x_target,y_target);
x_target = x_target(:);
y_target = y_target(:);
num_points_target = nx_target*ny_target;
Nt = size(z_array,2);
Nt_meas = round(T_meas/dt);
Nt_pred = round(T_pred/dt);
Nt_delay = round(T_delay/dt);
t_meas_start_ind = round(1:(Nt_meas*overlap):(Nt-Nt_meas+1));
t_pred_start_ind = t_meas_start_ind + Nt_meas + Nt_delay;
num_bursts = length(t_meas_start_ind);
t_pred = nan(num_bursts,Nt_pred,num_points_target);
z_target_pred = nan(num_bursts,Nt_pred,num_points_target);
x_target_vect = repmat(x_target(:)',[Nt_pred,1]);
y_target_vect = repmat(y_target(:)',[Nt_pred,1]);
n_points = length(x_array);

% loop over bursts, make least squares prediction
for i = 1:num_bursts
    % collect buoy array measurements for this burst
    t_target_ind = t_meas_start_ind(i):(t_meas_start_ind(i)+Nt_meas-1);
    x_array_vect = reshape(repmat(x_array(:)',[1,Nt_meas]),[1,Nt_meas*n_points]);
    y_array_vect = reshape(repmat(y_array(:)',[1,Nt_meas]),[1,Nt_meas*n_points]);
    z_array_vect = reshape(z_array(:,t_target_ind),[1,Nt_meas*n_points]);
    t_array_vect = reshape(repmat(Y.t(t_target_ind)',[n_points,1]),[1,Nt_meas*n_points]);
    % set time vector for prediction
    t_pred_ind = t_pred_start_ind(i):(t_pred_start_ind(i)+Nt_pred-1);
    t_pred(i,:,:) = repmat((t_pred_ind'-1)*dt,[1,1,num_points_target]);
    % least squares calculation
    [z_target_pred_vect,~] = leastSquaresWavePropagation_2D(z_array_vect,...
        t_array_vect,x_array_vect,y_array_vect,...
        t_pred(i,:,:),x_target_vect(:),y_target_vect(:),k,theta,reg_factor);
    z_target_pred(i,:,:) = reshape(z_target_pred_vect,[1,Nt_pred,num_points_target]);
end

% Populate "ground truth" comparison (from original simulation)
x_target_ind = nan(num_points_target,1);
y_target_ind = nan(num_points_target,1);
for i = 1:num_points_target
    x_target_ind(i) = find(Y.x == x_target(i));
    y_target_ind(i) = find(Y.y == y_target(i));
end
z_target_truth = nan(num_bursts,Nt_pred,num_points_target);
for i = 1:num_bursts
    t_pred_ind = t_pred_start_ind(i):(t_pred_start_ind(i)+Nt_pred-1);
    t_ind_compare = t_pred_ind(t_pred_ind<=Nt & t_pred_ind>=1);
    for j = 1:num_points_target
        z_target_truth(i,t_pred_ind<=Nt & t_pred_ind>=1,j) = squeeze(Y.Z(y_target_ind(j),x_target_ind(j),t_ind_compare));
    end
end