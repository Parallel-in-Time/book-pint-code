sigma=10;r=28;b=8/3;                  % Lorenz rhs and Jacobian
f=@(t,u) [sigma*(u(2)-u(1)) r*u(1)-u(2)-u(1)*u(3) u(1)*u(2)-b*u(3)];
jac=@(t,u) [[-sigma, sigma, 0    ]
            [r-u(3), -1   , -u(1)]
            [u(2)  , u(1) , -b   ]];
M=10;                                 % propagator of the coupled system
P=@(t0,t1,u0) CForwardEuler(f,jac,t0,t1,u0,M);
T=1; u0=[20;5;-5]; K=9; N=500;        % multiple shooting parameters
[~,uPred]=ForwardEuler(f,[0 T],u0,N); % compute initial guess using N
                                      % steps of Forward Euler
U=MultipleShooting(P,T,u0,N,K,uPred); % solve with multiple shooting
[t,u]=ForwardEuler(f,[0 T],u0,M*N);   % fine solution
TT=0:T/N:T;                           % coarse time mesh
for k=1:3
    plot3(u(:,1),u(:,2),u(:,3),'-b', ...         % plot fine solution
          U{k}(:,1),U{k}(:,2),U{k}(:,3),'.');    % and Multiple Shooting
    axis([-20 30 -30 40 -10 60]); view([-13,8]); % iterate
    xlabel('x'); ylabel('y'); zlabel('z');
    grid on
    pause
end