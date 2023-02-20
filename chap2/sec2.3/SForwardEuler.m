function u=SForwardEuler(f,t0,t1,u0,n);
[t,u]=ForwardEuler(f,[t0 t1],u0,n);
u=u(end,:);
