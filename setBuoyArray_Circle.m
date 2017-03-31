function [x_array,y_array,z_array] = setBuoyArray_Circle(x_center,y_center,theta,radius,dx,dy,Y)
 
n_radius = length(radius(:));
n_theta = length(theta(:));
radius_vec = repmat(radius(:),[1,n_theta]);
theta_vec = repmat(theta(:)',[n_radius,1]);
n_points = n_theta*n_radius;
x_array = x_center + round(radius_vec(:).*cos(theta_vec(:)*pi/180)/dx)*dx;
y_array = y_center + round(radius_vec(:).*sin(theta_vec(:)*pi/180)/dy)*dy;
x_array_ind = nan(n_points,1);
y_array_ind = nan(n_points,1);
Nt = size(Y.Z,3);
z_array = nan(n_points,Nt);
for i = 1:n_points
    x_array_ind(i) = find(Y.x == x_array(i));
    y_array_ind(i) = find(Y.y == y_array(i));
    z_array(i,:) = squeeze(Y.Z(y_array_ind(i),x_array_ind(i),:));
end