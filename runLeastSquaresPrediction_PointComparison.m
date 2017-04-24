function [z_target_pred,z_target_truth,t_pred] = runLeastSquaresPrediction_PointComparison(...
    x_target,y_target,z_target,dt,x_array,y_array,z_array,...
    k,theta,reg_factor,T_meas,T_pred,T_delay,overlap)


% reshape/initialize arrays
[Nt,numArray] = size(x_array);
Nt_meas = round(T_meas/dt);
Nt_pred = round(T_pred/dt);
Nt_delay = round(T_delay/dt);
t_meas_start_ind = round(1:(Nt_meas*overlap):(Nt-Nt_meas+1));
t_pred_start_ind = t_meas_start_ind + Nt_meas + Nt_delay;
num_bursts = length(t_meas_start_ind);
t_pred = nan(num_bursts,Nt_pred);
z_target_pred = nan(num_bursts,Nt_pred);
z_target_truth = nan(num_bursts,Nt_pred);

% loop over bursts, make least squares prediction
for i = 1:num_bursts
    % collect buoy array measurements for this burst
    t_target_ind = t_meas_start_ind(i):(t_meas_start_ind(i)+Nt_meas-1);
    x_array_vect = reshape(x_array(t_target_ind,:),[1,Nt_meas*numArray]);
    y_array_vect = reshape(y_array(t_target_ind,:),[1,Nt_meas*numArray]);
    z_array_vect = reshape(z_array(t_target_ind,:),[1,Nt_meas*numArray]);
    t_array_vect = reshape(repmat((t_target_ind'-1)*dt,[numArray,1]),[1,Nt_meas*numArray]);
    % set time vector for prediction
    t_pred_ind = t_pred_start_ind(i):(t_pred_start_ind(i)+Nt_pred-1);
    t_pred(i,:) = (t_pred_ind'-1)*dt;
    if max(t_pred_ind)<Nt
        % least squares calculation
        [z_target_pred_vect,~] = leastSquaresWavePropagation(z_array_vect,...
            t_array_vect,x_array_vect,y_array_vect,...
            t_pred(i,:),x_target(t_pred_ind,:),y_target(t_pred_ind,:),k,theta,reg_factor);
        z_target_pred(i,:) = reshape(z_target_pred_vect,[1,Nt_pred]);
        z_target_truth(i,:) = z_target(t_pred_ind);
    end
end
