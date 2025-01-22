c=1;                                           % wave speed
n=200;                                         % number of spatial points 
X=1;                                           % space domain (0,X)
dx=X/(n+1);                                    % spatial mesh size 
ex=ones(n,1);
DDx=1/dx^2*spdiags([ex -2*ex ex],[-1:1],n,n);  % second derivative in space
T=1;                                           % time domain (0,T) 
dtOnCFL=dx/c;                                  % dt on CFL, need to
m=ceil(T/dtOnCFL);                             % be just below for more  
dt=T/m;                                        % dt to be used
et=ones(m,1);
DDt=1/(c^2*dt^2)*spdiags([et -2*et et],[-2:0],m,m); 
                                               % second derivative in time
A=kron(DDt,speye(size(DDx)))-kron(spdiags(et,-1,m,m),DDx); 
                                               % all-at-once matrix
t=dt:dt:T;                                     % time and space mesh 
x=(dx:dx:1-dx)';
u0=exp(-200*(x-1/2).^2);                       % initial conditions
u0t=zeros(size(x));                            
f=zeros(size(A,1),1);                          % add initial conditions  
f(1:n)=u0t/dt+u0/dt^2+1/2*DDx*u0;              % on the right hand side
f(n+1:2*n)=-u0/dt^2;                           % using Taylor 
ue=A\f;                                        % exact solution 
Ue=reshape(ue,n,m); surf(t,x,Ue);              % reshape for plotting
xlabel('t');ylabel('x'); view(-90,90); colorbar
pause

nx=5;                                          % 5 red subdomains in space 
[Rr,Rb,mtr,mtb]=RedBlackSubdomains(n,m,nx);    % construct RAS matrices
u=rand(m*n,1)-1/2;                             % random initial guess
for jr=1:mtr                                   % run unmapped tent pitching
  for ir=1:nx                                  % red subdomain solves
    u=u+Rr{ir,jr}'*((Rr{ir,jr}*A*Rr{ir,jr}')\(Rr{ir,jr}*(f-A*u)));
  end;  
  U=reshape(ue-u,n,m);                         % to plot the error
  surf(t,x,U); xlabel('t');ylabel('x');        
  axis([0 T 0 X -1 1 -1 1]); colorbar; view(-90,90); pause
  if jr<=mtb                                   % also another black subdomain?
    for ib=1:nx-1                              % black subdomain solves
      u=u+Rb{ib,jr}'*((Rb{ib,jr}*A*Rb{ib,jr}')\(Rb{ib,jr}*(f-A*u)));
    end  
    U=reshape(ue-u,n,m);
    surf(t,x,U); xlabel('t');ylabel('x');
    axis([0 T 0 X -1 1 -1 1]); colorbar; view(-90,90); pause
  end;  
end;