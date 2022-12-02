function [u1,V]=CForwardEuler(f,jac,t0,t1,u0,n)
nDOF=length(u0);
u0=reshape([u0, eye(nDOF)],[],1);
rhs=@(t,u) [f(t,u(1:nDOF)); 
            reshape(jac(t,u(1:nDOF))*reshape(u(nDOF+1:end),nDOF,nDOF),[],1)];
[~,u]=ForwardEuler(rhs,[t0 t1],u0,n);
u=u(:,end);
u1=u(1:nDOF);
V=reshape(u(nDOF+1:end),nDOF,nDOF);

