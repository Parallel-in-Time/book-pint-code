lambda=1j-0.02;
T=6*pi;
M=5;
N=10;
u0=1;
[t,uIDC]=IDC(lambda,[0,T],u0,N,M,0);
uExact=exp(lambda*t);
plot(t,real(uExact),'o-'); hold on;
plot(t,real(uIDC))
for K=1:3
  [t,uIDC]=IDC(lambda,[0, T],u0,N,M,K);
  plot(t,real(uIDC))
end