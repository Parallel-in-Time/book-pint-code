% Settings for Matlab
close all; clear;
set( 0, 'defaultAxesTickLabelInterpreter', 'latex' );
set( 0, 'defaultLegendInterpreter',        'latex' );
set( 0, 'defaultTextInterpreter',          'latex' );
set( 0, 'defaultAxesFontSize', 12 );
format long 

% Newton parameters
x0 = 2;
nPas = 6;
global approx approxType h;
h = 1e-1;
vPas = 1:(nPas+1);

% Using exact derivative
approx = false;
approxType = 'NONE';
xk = zeros(1, nPas+1);
xk(1) = x0;
for k = 1:nPas
    xk(k+1) = xk(k) - f(xk(k))/fPrim(xk(k));
end
semilogy(vPas, abs(xk-pi), '-o', 'DisplayName', "f' exact");
hold on;

% Using approximate constant derivative
approx = true;
approxType = 'EVAL';
xk = zeros(1, nPas+1);
xk(1) = x0;
for k = 1:nPas
    xk(k+1) = xk(k) - f(xk(k))/fPrim(xk(k));
end
semilogy(vPas, abs(xk-pi), '-s', 'DisplayName', "f' const");

% Using finite difference approximation
approxType = 'FINDIFF';
xk = zeros(1, nPas+1);
xk(1) = x0;
for k = 1:nPas
    xk(k+1) = xk(k) - f(xk(k))/fPrim(xk(k));
end
semilogy(vPas, abs(xk-pi), '-^', 'DisplayName', "f' finDiff");

grid on; xlabel('$N_{iter}$'); ylabel('error'); legend()


function sol = f(x)
    sol = sin(x);
end

function sol = fPrim(x)
    global approx approxType h; 
    if approx
        if strcmp(approxType,'FINDIFF')
            sol = (f(x+h)-f(x-h))/(2*h);
        elseif strcmp(approxType, 'EVAL')
            sol = -0.90;
        end
    else
        sol = cos(x);
    end
end
