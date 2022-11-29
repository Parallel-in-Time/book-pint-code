f=@(x,t) zeros(size(x));                    % zero right hand side
u0=@(x) exp(-120*(x-0.3).^2);               % initial condition
g=@(t) zeros(size(t));                      % boundary condition
[u,x,t]=TransportFEUpwind(f,u0,g,0,1,40,0.5,20);
PlotTransport(u,x,t,u0);