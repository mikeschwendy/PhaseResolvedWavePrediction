function [z2,zc] = leastSquaresWavePropagation(z1,t1,x1,y1,t2,x2,y2,k,theta,reg_factor)
k = k(:);
theta = theta(:);
kx = k*cos(theta');
ky = k*sin(theta');
omega = sqrt(9.8*k)*ones(size(theta'));

kx = kx(:);
ky = ky(:);
omega = omega(:);
x1 = x1(:);
y1 = y1(:);
t1 = t1(:);
z1 = z1(:);
x2 = x2(:);
y2 = y2(:);
t2 = t2(:);

N_input_pts = length(z1);
if length(x1) ~= N_input_pts || length(y1) ~= N_input_pts || length(t1) ~= N_input_pts
    error('All input vectors must be equal length')
end

N_output_pts = length(t2);
if length(x2) ~= N_output_pts || length(y2) ~= N_output_pts
    error('All output vectors must be equal length')
end

% Note: Fitting complex exponentials (as in Connell et al. 2015) gives results
% that are off by a factor of 0.5.  It is not clear why.

% P1 = exp(1i.*(x1*kx'+y1*ky'-t1*omega'));
% A = ridge(z1,P1,reg_factor,0);
% zc = real(P1*A(2:end));
% P2 = exp(1i.*(x2*kx'+y2*ky'-t2*omega'));
% z2 = real(P2*A(2:end));

% However, fitting cosines and sines works.  See test script below for
% proof.
P1 = [cos(x1*kx'+y1*ky'-t1*omega'),sin(x1*kx'+y1*ky'-t1*omega')];
A = ridge(z1,P1,reg_factor,0);
zc = P1*A(2:end);
P2 = [cos(x2*kx'+y2*ky'-t2*omega'),sin(x2*kx'+y2*ky'-t2*omega')];
z2 = P2*A(2:end);

%% Test script
% Script shows the difference between fitting complex exponentials vs.
% sines and cosines
if false
    x = linspace(0,10,301);
    k = linspace(1,10,10);
    z = 5*sin(x*3)+2*cos(x*4);
    P1 = exp(1i.*(x'*k));
    P2 = [cos(x'*k),sin(x'*k)];
    A1 = P1\z';
    A2 = P2\z';
    zc1 = real(P1*A1);
    zc2 = P2*A2;
    
    figure(1); clf(1);
    subplot(2,1,1)
    plot(x,z)
    hold on
    plot(x,zc1,'--')
    plot(x,zc2,'-.')
    hold off
    legend('Original','Exponential Fit','Sin/cos Fit')
    subplot(2,1,2)
    plot(k,real(A1))
    hold on
    plot(k,imag(A1))
    plot(k,A2(1:10))
    plot(k,A2(11:20))
    hold off
    legend('Real A1','Imag A1','A2cos()','A2sin()')
end
