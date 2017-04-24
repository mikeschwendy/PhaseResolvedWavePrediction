function [z_target_pred,t_pred] = runLeastSquaresPrediction_SurfaceReconstruction(...
    x_target,y_target,x_array,y_array,z_array,dt,...
    k,theta,reg_factor,T_meas,T_pred,T_delay,overlap)


% reshape/initialize arrays
x_target = x_target(:);
y_target = y_target(:);
num_points_target = length(x_target);
[Nt_tot,num_points_array] = size(x_array);
Nt_meas = round(T_meas/dt);
Nt_pred = round(T_pred/dt);
Nt_delay = round(T_delay/dt);
t_meas_start_ind = round(1:(Nt_meas*overlap):(Nt_tot-Nt_meas+1));
t_pred_start_ind = t_meas_start_ind + Nt_meas + Nt_delay;
num_bursts = length(t_meas_start_ind);
t_pred = nan(num_bursts,Nt_pred,num_points_target);
z_target_pred = nan(num_bursts,Nt_pred,num_points_target);
x_target_vect = repmat(x_target(:)',[Nt_pred,1]);
y_target_vect = repmat(y_target(:)',[Nt_pred,1]);

% loop over bursts, make least squares prediction
for i = 1:num_bursts
    % collect buoy array measurements for this burst
    t_target_ind = t_meas_start_ind(i):(t_meas_start_ind(i)+Nt_meas-1);
    x_array_vect = reshape(x_array(t_target_ind,:),[1,Nt_meas*num_points_array]);
    y_array_vect = reshape(y_array(t_target_ind,:),[1,Nt_meas*num_points_array]);
    z_array_vect = reshape(z_array(t_target_ind,:),[1,Nt_meas*num_points_array]);
    t_array_vect = reshape(repmat((t_target_ind'-1)*dt,[num_points_array,1]),[1,Nt_meas*num_points_array]);
    % set time vector for prediction
    t_pred_ind = t_pred_start_ind(i):(t_pred_start_ind(i)+Nt_pred-1);
    t_pred(i,:,:) = repmat((t_pred_ind'-1)*dt,[1,1,num_points_target]);
    % least squares calculation
    [z_target_pred_vect,~] = leastSquaresWavePropagation(z_array_vect,...
        t_array_vect,x_array_vect,y_array_vect,...
        t_pred(i,:,:),x_target_vect(:),y_target_vect(:),k,theta,reg_factor);
    z_target_pred(i,:,:) = reshape(z_target_pred_vect,[1,Nt_pred,num_points_target]);
end
