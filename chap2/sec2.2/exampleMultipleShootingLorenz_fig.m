% Settings for Matlab
close all; clear;
set(0, 'defaultLegendInterpreter', 'latex');
format long

sigma=10;r=28;b=8/3;                        % Lorenz rhs and Jacobian
f=@(t,u) [sigma*(u(2)-u(1)) 
          r*u(1)-u(2)-u(1)*u(3) 
          u(1)*u(2)-b*u(3)];
jac=@(t,u) [[-sigma, sigma, 0    ]
            [r-u(3), -1   , -u(1)]
            [u(2)  , u(1) , -b   ]];
M=10;
propagator=@(t0,t1,u0) CForwardEuler(f,jac,t0,t1,u0,M);
u0=[20;5;-5]; K=9; N=500;

for T=[0.5, 1, 1.5, 1.78, 2]

    [tFine,uFine]=ForwardEuler(f,[0 T],u0,M*N);
    uRef=uFine(:,1:M:end);
    [tPred,uPred]=ForwardEuler(f,[0 T],u0,N);
    [t,U]=MultipleShooting(propagator,T,u0,N,K,uPred);

    err = squeeze(vecnorm(vecnorm(U-uRef,2),Inf));

    figure(1); 
    semilogy(err, DisplayName=sprintf('$T=%1.2f$', T)); 
    hold on; xlabel('k'); ylabel('error');
    xlim([1 10]); ylim([1e-14 1e7]);
    legend(Location='east'); set(gca,'fontsize', 12);
    
    if (T == 1 || T == 2)
        errTime = squeeze(vecnorm(U-uRef,2));
        errTime(errTime==0) = 1e-16;
        figure(T+1); 
        for k=1:7
            semilogy(t, errTime(:, k), DisplayName=sprintf('$k=%d$', k-1));
            hold on; xlabel('time'); ylabel('error');
            xlim([-T/20, T]); ylim([1e-15, 1e4]);
            legend(Location='southeast'); set(gca,'fontsize', 12);
        end
    end   
end

% Save Figures into files
figure(1); 
print -depsc ../../../book-pint/Figures/MultipleShootingLorenzError.eps;
figure(2); 
print -depsc ../../../book-pint/Figures/MultipleShootingLorenzError_T=1.eps;
figure(3); 
print -depsc ../../../book-pint/Figures/MultipleShootingLorenzError_T=2.eps;