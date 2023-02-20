sigma=10;r=28;b=8/3;                                % Lorenz rhs
f=@(t,x) [sigma*(x(2)-x(1)) r*x(1)-x(2)-x(1)*x(3) x(1)*x(2)-b*x(3)];
MF=10; MG=1;                                        % F and G time steps
F=@(t0,t1,u0) SForwardEuler(f,t0,t1,u0,MF);         % fine solver F
G=@(t0,t1,u0) SForwardEuler(f,t0,t1,u0,MG);         % coarse solver G
K=20; u0=[20;5;-5]; T=5; N=500;                          
U=Parareal(F,G,T,u0,N,K);                                
[t,u]=ForwardEuler(f,[0 T],u0,MF*N);                % fine solution
TT=0:T/N:T;
for k=1:K                                           % plot fine
  plot3(u(:,1),u(:,2),u(:,3),'-b'...                % solution and
    ,U{k}(:,1),U{k}(:,2),U{k}(:,3),'.');            % parareal iterate
  axis([-20 30 -30 40 -10 60]); view([-13,8]);
  xlabel('x'); ylabel('y'); zlabel('z');
  set(gca,'Fontsize',18)
  grid on 
  pause
end

for k=1:K+1                             
  err(k)=max(max(abs(u(1:10:end,:)-U{k})))
end      

semilogy(0:K,err,'-')
set(gca,'Fontsize',12)
xlabel('k')          
ylabel('error')
legend('parareal error')
print -depsc PararealLorenzError.eps
