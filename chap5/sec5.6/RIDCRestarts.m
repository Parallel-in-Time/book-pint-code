function [t,u]=RIDCRestarts(lambda,tspan,u0,N,M,K,R)
% RIDCRESTARTS Solves the Dahlquist problem using RIDC with restarts.
% [t,u]=RIDCRestarts(lambda,tspan,u0,N,M,K,R) solves the Dahlquist
    test equation using RIDC with R restarts.

NR=N/R; t0=tspan(1); dtR=(tspan(2)-tspan(1))/R;
uR0=u0; t=tspan(1); u=u0;
for r=1:R
  tspanR=[t0+(r-1)*dtR,t0+r*dtR];
  [tR,uR]=RIDC(lambda,tspanR,uR0,NR,M,K);
  t=[t,tR(2:end)];
  u=[u,uR(2:end)];
  uR0=u(end);
end
