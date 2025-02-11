sigma=10;r=28;b=8/3;
f=@(t,x) [sigma*(x(2)-x(1)); r*x(1)-x(2)-x(1)*x(3); x(1)*x(2)-b*x(3)];
T=30;N=30000;dt=T/N;
x0=[20;5;-5];
[t,xS]=MirankerLinigerS(f,[0 T],x0,N);
plot3(xS(1,:),xS(2,:),xS(3,:),'-b');
