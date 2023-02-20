f=@(x,t) x.^4.*(1-x).^4+10*sin(8*t);
%f=@(x,t) 100*exp(-100*(x-0.5).^2)*sin(t);
%f=@(x,t) sin(pi*x).*100*sin(t);
%f=@(x,t) 0.2+sin(2*pi*x).*sin(0.5*t);

%f=@(x,t) 100*exp(-50*((x-0.1).^2)-5*(t-1).^2)-15*exp(-30*((x-0.6).^2)-5*(t-1).^2)+100*exp(-50*((x-0.9).^2)-5*(t-3).^2)-15*exp(-30*((x-0.4).^2)-5*(t-3).^2)+100*exp(-50*((x-0.1).^2)-5*(t-5).^2)-15*exp(-30*((x-0.6).^2)-5*(t-5).^2)+100*exp(-50*((x-0.9).^2)-5*(t-7).^2)-15*exp(-30*((x-0.4).^2)-5*(t-7).^2);

T=8; N=16; K=16; J=10; 

u0=zeros(J+1,1);
MG=1;                                             % MF no of fine steps
gl=zeros(MG+1,1);gr=zeros(MG+1,1);
G=@(t0,t1,u0) SHeatEquationBE(f,[t0 t1],[0 1],u0,gl,gr);  % G coarse solver
MF=10;                                            % MF no of fine steps
gl=zeros(MF+1,1);gr=zeros(MF+1,1);
F=@(t0,t1,u0) SHeatEquationBE(f,[t0 t1],[0 1],u0,gl,gr);  % F fine solver
%F=@(t0,t1,u0) SHeatEquationFE(f,[t0 t1],[0 1],u0,gl,gr);  % F fine solver

U=Parareal(F,G,T,u0,N,K);                         
glf=zeros(MF*N+1,1); grf=zeros(MF*N+1,1); 
u=HeatEquationBE(f,[0 T],[0 1],u0,glf,grf);       % fine solution
%u=HeatEquationFE(f,[0 T],[0 1],u0,glf,grf);       % fine solution
dt=T/(MF*N);dT=T/N; dx=1/J;                       % mesh parameters
t=(0:dt:T)'; TT=(0:dT:T)'; x=0:dx:1;
Uf=u(1:MF:end,:);
subplot(1,2,1)
mesh(x,t,u,'linewidth',2); xlabel('x'); ylabel('t'); title('solution')
axis([0 1 0 T -1 1])
set(gca,'fontsize',24)
pause
%print('-depsc','PararealHeatSolution.eps')
for k=1:K
  up=[];
  for n=1:N
    up((n-1)*MF+1:n*MF+1,:)=HeatEquationBE(f,[(n-1)*dT n*dT],[0 1],U{k}(n,:),gl,gr);
    %up((n-1)*MF+1:n*MF+1,:)=HeatEquationFE(f,[(n-1)*dT n*dT],[0 1],U{k}(n,:),gl,gr);
  end;  
  up(N*MF+1,:)=U{k}(end,:);
  subplot(1,2,1)
  mesh(x,t,up,'linewidth',2); xlabel('x'); ylabel('t');title(['iteration ' num2str(k)]) 
  set(gca,'Fontsize',24)
  axis([0 1 0 T -1 1])
  %print('-depsc',['PararealHeatIter' num2str(k) '.eps'])
  subplot(1,2,2)
  mesh(x,t,u-up); xlabel('x'); ylabel('t');title('error')
  set(gca,'Fontsize',24)
  %print('-depsc',['PararealHeatErrorIter' num2str(k) '.eps'])
  pause
  err(k)=max(max(abs(u-up)));
end  

semilogy(1:K,err,'-','linewidth',2)
set(gca,'Fontsize',24)
xlabel('k')          
ylabel('error')
legend('parareal error')
%print -depsc PararealHeatErrorCurve.eps



