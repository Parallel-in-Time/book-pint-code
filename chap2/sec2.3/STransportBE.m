function u=STransportBE(f,a,tspan,xspan,u0,N);
u=TransportBE(f,a,tspan,xspan,u0,N);
u=u(end,:);