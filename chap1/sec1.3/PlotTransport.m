function PlotTransport(u,x,t,u0);
% PLOTTRANSPORT plots an numerical and exact transport solution
%   PlotTransport(u,x,t,u0); plots a numerical solution in u on the
%   space grid in x and time grid in t (colums of u are space),
%   representing an approximate solution for an advection equation,
%   and also the exact solution from the initial condition is plotted
%   if given by the initial condition function u0.

ma=max(max(u));
mi=min(min(u));
a=x(1);
b=x(length(x));
xx=(a:(b-a)/500:b)';
for n=1:size(u,2)
  if nargin==4
    plot(x,u(:,n),'o',xx,feval(u0,xx-t(n)),'-','LineWidth',2,'Markersize',10);
  else      
    plot(x,u(:,n),'o')
  end  
  xlabel('x');ylabel('u');
  t(n)
  pause
end;