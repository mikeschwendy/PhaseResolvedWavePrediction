function [z2,zc] = leastSquaresWavePropagation(z1,t1,x1,t2,x2,f)

%f_ls = [-fSim(11:2:100)+1e-3,fSim(10:2:100)];
k = -sign(f).*(2*pi*f).^2/(9.8);
Nf = length(f);

Nx1 = length(x1);
Nt1 = length(t1);
k1 = repmat(k(:)',[Nx1*Nt1,1]);
f1 = repmat(f(:)',[Nx1*Nt1,1]);
t1 = repmat(t1(:),[1,Nf]);
x1 = repmat(x1(:),[1,Nf]);

Nx2 = length(x2);
Nt2 = length(t2);
k2 = repmat(k(:)',[Nx2*Nt2,1]);
f2 = repmat(f(:)',[Nx2*Nt2,1]);
t2 = repmat(t2(:),[1,Nf]);
x2 = repmat(x2(:),[1,Nf]);

P1 = exp(1i.*(x1.*k1-t1.*f1*2*pi));
A = P1\z1(:); 

zc = reshape(real(P1*A),[Nx1,Nt1]);

P2 = exp(1i.*(x2.*k2-t2.*f2*2*pi));
z2 = reshape(real(P2*A),[Nx2,Nt2]);



