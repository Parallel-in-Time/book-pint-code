l=7;                                             % with the choice l=3 and
J=2^l-1;                                         % N=(J+1)*2^(3*l)-1 the 
N=J;                                             % space coarsening condition 
e=ones(J,1);                                     % fails, try it 
dx=1/(J+1); x=(0:dx:1)';
A=1/dx^2*spdiags([e -2*e e],[-1 0 1],J,J);
T=5; dt=T/N; t=(0:dt:T);
u0=@(x) zeros(size(x));
gl=@(t) zeros(size(t)); gr=@(t) zeros(size(t));
f=@(x,t) x.^4.*(1-x).^4+10*sin(8*t);
for n=1:N                                        % compute rhs b for reuse
  b(:,n)=dt*feval(f,x(2:end-1),t(n+1));          % source function
  b(1,n)=b(1,n)+dt/dx^2*gl(t(n+1));              % boundary conditions
  b(end,n)=b(end,n)+dt/dx^2*gr(t(n+1));
end;

u(:,1)=u0(x); u(1,1:N+1)=gl(t); u(end,1:N+1)=gr(t);
G=speye(size(A))-dt*A;
for n=1:N                                        % compute exact solution
  u(2:end-1,n+1)=G\(u(2:end-1,n)+b(:,n));        % exact BE
end
uBE=u;

Jc=(J+1)/2-1;                                    % coarse grid size
P=sparse(J,Jc);                                  % prolongation by interpolation
for j=1:Jc
  P(2*j,j)=1; P(2*j-1,j)=0.5; P(2*j+1,j)=0.5;
end;
R=0.5*P';                                        % restriction
Ac=R*A*P;                                        % coarse matrix by Galerkin 
Nc=(N+1)/2-1;
Pt=sparse(N,Nc);                                 % prolongation by interpolation
for j=1:Nc
  Pt(2*j,j)=1; Pt(2*j-1,j)=0.5; Pt(2*j+1,j)=0.5;
end;
Rt=0.5*Pt';                                      % restriction
Pt(end,end)=1;                                   % no final zero bc in time
Gc=speye(size(Ac))-2*dt*Ac;                      % coarsening in time also

rng('default')                                   % random initial guess
u(2:end-1,2:end)=rand(J,N);                      % with correct ic and bc
errcxt(1)=max(max(abs(uBE-u)));
for k=1:10
  for j=1:nu                                     % use nu block Jacobi steps
    uo=u;                                        % block Jacobi is parallel
    for n=1:N
      u(2:end-1,n+1)=(1-al)*uo(2:end-1,n+1)+al*(G\(uo(2:end-1,n)+b(:,n)));
    end;
  end
  mesh(x,t,uBE'-u'); xlabel('x'); ylabel('t');
  title('Error after presmoothing');
  pause
  for n=1:N                                      % compute residual
    r(:,n)=u(2:end-1,n)+b(:,n)-G*u(2:end-1,n+1);
  end;
  rc=R*r;                                        % restrict residual in space
  rc=rc*Rt';                                     % restrict residual in time 
  uc=zeros(Jc+2,1);                              % zero ic for correction
  for n=1:Nc                                     % coarse correction
    uc(2:end-1,n+1)=Gc\(uc(2:end-1,n)+rc(:,n));  % exact BE using 2*dt
  end                                            % extend in space and time
  u(2:end-1,2:end)=u(2:end-1,2:end)+P*uc(2:end-1,2:end)*Pt';
  mesh(x,t,uBE'-u'); xlabel('x'); ylabel('t');
  title('Error after coarse correction');
  pause
  for j=1:nu                                     % use nu block Jacobi steps
    uo=u;                               
    for n=1:N
      u(2:end-1,n+1)=(1-al)*uo(2:end-1,n+1)+al*(G\(uo(2:end-1,n)+b(:,n)));
    end;
  end
  errcxt(k+1)=max(max(abs(uBE-u)));
  mesh(x,t,uBE'-u'); xlabel('x'); ylabel('t');
  title(['Error after 2 grid iteration k=' num2str(k)]);
  pause
end


