f=@(x,t) x.^4.*(1-x).^4+10*sin(8*t);              % heat source function
T=8; N=16; K=16; J=10;                            % parareal parameters
u0=zeros(J+1,1);                                  % initial condition
MG=1;                                             % MG no of coarse steps
gl=zeros(MG+1,1);gr=zeros(MG+1,1);                % G boundary conditions
G=@(t0,t1,u0) SHeatEquationBE(f,[t0 t1],[0 1],u0,gl,gr);  % G coarse solver
MF=10;                                            % MF no of fine steps
gl=zeros(MF+1,1);gr=zeros(MF+1,1);                % F boundary conditions
F=@(t0,t1,u0) SHeatEquationBE(f,[t0 t1],[0 1],u0,gl,gr);  % F fine solver
U=Parareal(F,G,T,u0,N,K);                         
glf=zeros(MF*N+1,1); grf=zeros(MF*N+1,1);         
u=HeatEquationBE(f,[0 T],[0 1],u0,glf,grf);       % fine solution
dt=T/(MF*N);dT=T/N; dx=1/J;                       % mesh parameters
t=(0:dt:T)'; TT=(0:dT:T)'; x=0:dx:1;              % for plotting
for k=1:K
  up=[];
  for n=1:N                                       % reconstruct fine
    up((n-1)*MF+1:n*MF+1,:)=HeatEquationBE(f,...  % solution from 
        [(n-1)*dT n*dT],[0 1],U{k}(n,:),gl,gr);   % parareal for plotting
  end;  
  up(N*MF+1,:)=U{k}(end,:);
  mesh(x,t,up); xlabel('x'); ylabel('t');         % plot parareal approx.
  axis([0 1 0 T -1 1])
  pause
  mesh(x,t,u-up); xlabel('x'); ylabel('t');       % plot parareal error
  err(k)=max(max(abs(u-up)));
  pause
end