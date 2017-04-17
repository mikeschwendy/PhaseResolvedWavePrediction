function [x_array,y_array,z_array] = setBuoyArray_Box(x_center,y_center,box_side_length,num_buoys,dir_principal,dx,dy,Y)

if num_buoys == 4
    x_array = x_center + box_side_length/2*[-1 0 1 0];
    y_array = y_center + box_side_length/2*[0 1 0 -1];
elseif num_buoys == 8
    x_array = x_center + box_side_length/2*[-1 0 1 0 0.5 0 -0.5 0];
    y_array = y_center + box_side_length/2*[0 1 0 -1 0 0.5 0 -0.5];
end

R = [cos(dir_principal*pi/180), -sin(dir_principal*pi/180); sin(dir_principal*pi/180), cos(dir_principal*pi/180)];
xy_array_rot = R*[(x_array-x_center);(y_array-y_center)];
x_array = round((x_center+xy_array_rot(1,:))/dx)*dx;
y_array = round((y_center+xy_array_rot(2,:))/dy)*dy;

x_array_ind = nan(num_buoys,1);
y_array_ind = nan(num_buoys,1);
Nt = size(Y.Z,3);
z_array = nan(num_buoys,Nt);
for i = 1:num_buoys
    x_array_ind(i) = find(Y.x == x_array(i));
    y_array_ind(i) = find(Y.y == y_array(i));
    z_array(i,:) = squeeze(Y.Z(y_array_ind(i),x_array_ind(i),:));
end

