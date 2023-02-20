f=@(x,t) 0;
a=1;
T=4; N=16; K=16;                                 
J=20; 
dx=1/J; x=0:dx:1;                                 % mesh parameters
u0=sin(2*pi*x);
MG=1;                                             % MF no of fine steps
G=@(t0,t1,u0) STransportBE(f,a,[t0 t1],[0 1],u0,MG); % G coarse solver
MF=20;                                            % MF no of fine steps
F=@(t0,t1,u0) STransportBE(f,a,[t0 t1],[0 1],u0,MF); % F fine solver
U=Parareal(F,G,T,u0,N,K);                         
u=TransportBE(f,a,[0 T],[0 1],u0,N*MF);             % fine solution
dt=T/(MF*N);dT=T/N;
t=(0:dt:T)'; TT=(0:dT:T)';
Uf=u(1:MF:end,:);
subplot(1,2,1);
mesh(x,t,u,'linewidth',2); xlabel('x'); ylabel('t');
axis([0 1 0 T -1 1])
set(gca,'fontsize',24)
pause
% print('-depsc','PararealTransportSolution.eps')
for k=1:K
  up=[];
  for n=1:N
    up((n-1)*MF+1:n*MF+1,:)=TransportBE(f,a,[(n-1)*dT n*dT],[0 1],U{k}(n,:),MF);
  end;  
  up(N*MF+1,:)=U{k}(end,:);
  subplot(1,2,1);
  mesh(x,t,up,'linewidth',2); xlabel('x'); ylabel('t');
  set(gca,'Fontsize',24)
  axis([0 1 0 T -1 1])
%  pause
%  print('-depsc',['PararealTransportIter' num2str(k) '.eps'])
  subplot(1,2,2);
  mesh(x,t,u-up,'linewidth',2); xlabel('x'); ylabel('t');
  set(gca,'Fontsize',24)
%  print('-depsc',['PararealTransportErrorIter' num2str(k) '.eps'])
  pause
  err(k)=max(max(abs(u-up)));
end  

semilogy(1:K,err,'-','linewidth',2)
set(gca,'Fontsize',24)
axis([0 16 1e-15 10])
xlabel('k')          
ylabel('error')
legend('parareal error a=0.25','parareal error a=0.5','parareal error a=1')
%print -depsc PararealTransportErrorCurve.eps
