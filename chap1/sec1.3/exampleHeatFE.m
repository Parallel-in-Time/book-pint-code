%N=1980;                                 % N number of time steps
N=2000;                                 % N number of time steps
J=100;                                  % J number of spatial steps 
e=ones(J-1,1);                          % N=2000 limit for J=100, T=1/10 
dx=1/J;                                 % finite difference Laplacian 
A=spdiags([e -2*e e],[-1 0 1],J-1,J-1)/dx^2;  
T=1/10;
dt=T/N;
u=20*e;                                 % initial temperature
for n=1:N                              
  u(:,n+1)=u(:,n)+dt*A*u(:,n);          % Forward Euler Step
  plot(0:dx:1,[0;u(:,n+1);0]);
  axis([0 1 0 21]); xlabel('x'); ylabel('u');
  pause
end;

