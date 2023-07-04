N=10; e=ones(N,1);                              % problem parameters
T=1; dt=T/N; t=0:dt:T;
la=-1; u0=1;
A=spdiags([-e/dt e/dt-la],[-1 0],N,N);          % time stepping matrix
f=zeros(N,1); f(1)=u0/dt;
ue=A\f;                                         % exact solution 
al=0.1; At=A; At(1,end)=al*At(2,1);             % make alpha-cyclic
D=spdiags(al.^((0:N-1)'/N),0,N,N);              % diagonal scaling
L=spdiags(conj(fft([0;1;zeros(N-2,1)])),0,N,N); % eigenvalues
u=zeros(N,1); up=u;                             % initial guess
K=5;
for k=1:K+1
  plot(t,[u0;real(ue)],'-',t,[u0;real(u)],'o',t,[u0;real(up)],'*');
  xlabel('t'); legend('solution','direct solve','ParaDiag solve');
  err(k)=(norm(ue-u));errp(k)=(norm(ue-up));    % compute errors
  pause
  u=u+At\(f-A*u);                               % u direct and up ParaDiag  
  up=up+dt*al^(-1/N)*inv(D)*fft((al^(-1/N)*(1-la*dt)*eye(N)-L)\ifft(D*(f-A*up)));
end;