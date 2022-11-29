sigma=10;r=28;b=8/3;
f=@(t,x) [sigma*(x(2)-x(1)); r*x(1)-x(2)-x(1)*x(3); x(1)*x(2)-b*x(3)];
T=30;N=30000;dt=T/N;
x=[20;5;-5];
for i=1:N
  x(:,i+1)=x(:,i)+dt*f(i*dt,x(:,i));             % Forward Euler step
  if mod(i,100)==0                               % plot only every 100th
    plot3(x(1,:),x(2,:),x(3,:),'-b');            % for animation speed
    axis([-20 30 -30 40 -10 60]); view([-13,8]);
    xlabel('x'); ylabel('y'); zlabel('z');
    grid on
    pause
  end  
end
hold on
xf=sqrt(b*(r-1)); yf=sqrt(b*(r-1)); zf=r-1;
plot3(xf,yf,zf,'or'); plot3(-xf,-yf,zf,'or');    % plot fixed points
plot3(0,0,0,'or'); 
hold off