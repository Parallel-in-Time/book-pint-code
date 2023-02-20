f=@(x,t) 0;
a=1;
T=4; N=16; K=16; J=20;                               % parareal parameters
dx=1/J; x=0:dx:1;                                    % spatial mesh
u0=sin(2*pi*x);                                      % initial condition
MG=1;                                                % MG no of coarse steps
G=@(t0,t1,u0) STransportBE(f,a,[t0 t1],[0 1],u0,MG); % G coarse solver
MF=20;                                               % MF no of fine steps
F=@(t0,t1,u0) STransportBE(f,a,[t0 t1],[0 1],u0,MF); % F fine solver
U=Parareal(F,G,T,u0,N,K);                          
u=TransportBE(f,a,[0 T],[0 1],u0,N*MF);              % fine solution
dt=T/(MF*N);dT=T/N; t=(0:dt:T)'; TT=(0:dT:T)';
for k=1:K
  up=[];
  for n=1:N                                          % reconstruct fine 
    up((n-1)*MF+1:n*MF+1,:)=TransportBE(f,a,...      % solution from 
             [(n-1)*dT n*dT],[0 1],U{k}(n,:),MF);    % parareal for plotting
  end;  
  up(N*MF+1,:)=U{k}(end,:);
  mesh(x,t,up); xlabel('x'); ylabel('t');            % plot parareal approx.
  axis([0 1 0 T -1 1])
  pause
  mesh(x,t,u-up); xlabel('x'); ylabel('t');          % plot parareal error
  err(k)=max(max(abs(u-up)));
  pause
end  