N=5000;                                 % 500 limit for M=50 and T=10
M=50;
e=ones(M-1,1);
L=1;
dx=L/M;
T=10;
dt=T/N;
x=0:dx:L;
for j=1:length(x)                                   % guitar initial condition
  if x(j)<L/2 u(j,1)=x(j); else u(j,1)=L-x(j); end;
  u(j,2)=u(j,1)+dt*0;                               % ut0=0 here
end;
%u(:,1)=0.6*exp(-100*(x-0.3).^2);
%u(:,1)=ones(size(u(:,1)));
%u(end,1:N+1)=1
%u(:,2)=u(:,1);   % other initial condition
for n=2:N 
  plot(x,u(:,n));
  axis([0 1 -0.6 0.6]);
  xlabel('x');ylabel('u');
  (n-1)*dt
  pause
  u(2:end-1,n+1)=2*u(2:end-1,n)-u(2:end-1,n-1)...
      +dt^2/dx^2*(u(1:end-2,n)-2*u(2:end-1,n)+u(3:end,n));
end;