function u=SDahlquistBE(la,t0,t1,u0,n);
[t,u]=DahlquistBE(la,[t0 t1],u0,n);
u=u(end);