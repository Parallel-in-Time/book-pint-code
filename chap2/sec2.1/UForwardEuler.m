function u=UForwardEuler(f,t0,t1,u0,n)
[~,u]=ForwardEuler(f,[t0 t1],u0,n);