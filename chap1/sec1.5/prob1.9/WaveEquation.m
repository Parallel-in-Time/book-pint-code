function [t, u] = WaveEquation(u0, ut0, c, L, tSpan, nStep)
    nDOF = length(u0);

    u = zeros(nDOF, nStep+1);
    t = linspace(0, tSpan, nStep+1);

    dx = L/(nDOF-1);
    dt = tSpan/nStep;

    u(:, 1) = u0;
    u(:, 2) = u0 + dt*ut0;

    coeff = c*dt^2/dx^2;
    for n =2:nStep
        u(2:end-1, n+1) = 2*u(2:end-1, n) - u(2:end-1, n-1) + ...
            coeff*(u(3:end, n) - 2*u(2:end-1, n) + u(1:end-2, n));
    end
end