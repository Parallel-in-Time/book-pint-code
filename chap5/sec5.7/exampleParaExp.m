T=0.4;                                  % time interval length 
g=@(x,t) 10*ones(size(x));              % heat source function
p=4;                                    % number of time intervals
N=10;                                   % N number of time steps per interval
J=10;                                   % J number of spatial steps 
e=ones(J-1,1);                            
dx=1/J; x=(dx:dx:1-dx)';                % finite difference Laplacian 
A=spdiags([e -2*e e],[-1 0 1],J-1,J-1)/dx^2;  
dt=T/(p*N); dT=T/p;                     % fine and coarse time step
u0=ones(J-1,1);                         % initial temperature
u(:,1)=u0;                              % compute reference solution
for n=1:N*p
  u(:,n+1)=(speye(J-1)-dt*A)\(u(:,n)+dt*g(x,dt*(n+1)));
end;
mesh(x,0:dt:T,u'); xlabel('x'); ylabel('t');
pause; hold on                          % keep to overlay images

for j=1:p                               % compute v solutions
  v{j}(:,1)=zeros(J-1,1); 
  for n=1:N                            
    v{j}(:,n+1)=(speye(J-1)-dt*A)\(v{j}(:,n)+dt*g(x,dt*((j-1)*N+n+1)));
  end;
  mesh(x,(j-1)*N*dt:dt:j*N*dt,v{j}'); pause
end;  

w{1}=u0;                                % compute w solutions
for j = 1:p
  for n=j:p
    w{j}(:,n+1)=expm(dT*A)*w{j}(:,n);   % should use efficient matrix exp
  end;  
  surf(x,0:dT:T,w{j}'); pause
  w{j+1}(:,j+1)=v{j}(:,N+1);
end;

UPE=w{1};                               % sum ParaExp solution 
surf(x,0:dT:T,UPE'); pause;
for j=2:p
  UPE=UPE+w{j};
  surf(x,0:dT:T,UPE'); pause
end;
UPE(:,p+1)=UPE(:,p+1)+v{p}(:,N+1);
surf(x,0:dT:T,UPE');
hold off