function u=SHeatEquationBE(f,tspan,xspan,u0,gl,gr);
u=HeatEquationBE(f,tspan,xspan,u0,gl,gr);
u=u(end,:);