%% Conductivity according to the Effective Medium Theory
% Depending on eta = n/nrms (mean number of impurities over root mean square of the
% distribution of impurities), mu (mobility), and B field.

function [sigmaxxemt, sigmaxyemt] = sigmaemt(eta, mu, B, nshift)

%% Conductivities
sigmaxx = @(nn, BB) sigmasimple(nn, mu, nshift) ./ (1 + (mu * BB) .^ 2);
sigmaxy = @(nn, BB) -sign(nn-nshift).*mu * BB .* sigmaxx(nn, BB);

denom = @(z, BB, y) (sigmaxx(z , BB) + y(1)) .^ 2 + (sigmaxy(z , BB) - y(2)) .^ 2;
%% Calculation Body
%For sc iteration
numiter1 = @(z, BB, y) (sigmaxx(z, BB) .^ 2 + (sigmaxy(z, BB) - y(2)) .^ 2) ./ denom (z, BB, y);
numiter2 = @(z, BB, y) sigmaxy(z, BB) ./ denom (z, BB, y);

%Integration limit
f = 5.3; %Optimized

%Gaussian distribution
p = @(z, z0) exp(-1 * ((z - z0) .^ 2) / 2) / ((sqrt(2 * pi)));

numiterp1 = @(z, z0, BB, y) p(z, z0) .* numiter1(z, BB, y);
numiterp2 = @(z, z0, BB, y) p(z, z0) .* numiter2(z, BB, y);
denomp = @(z, z0, BB, y) p(z, z0) ./ denom(z, BB, y); 

fnum1 = @(z0, BB, y) integral(@(z) numiterp1(z, z0, BB, y), z0 - f, z0 + f);
fnum2 = @(z0, BB, y) integral(@(z) numiterp2(z, z0, BB, y), z0 - f, z0 + f);
fdenom = @(z0, BB, y) integral(@(z) denomp(z, z0, BB, y), z0 - f, z0 + f);

fiter1 = @(z0, BB, y) fnum1(z0, BB, y) / fdenom(z0, BB, y);
fiter2 = @(z0, BB, y) fnum2(z0, BB, y) / fdenom(z0, BB, y);

%First guess based on sigma_min
nstar = 1 / sqrt(2);
sigmamin = sigmasimple(nstar, mu, nshift);
fsigmacrudexx = @(nn, BB) (sigmasimple(nn, mu, nshift) .* (abs(nn - nshift) > nstar) + sigmamin * (abs(nn - nshift) <= nstar)) ./ (1 + (mu * BB) .^ 2);

%% Solver
    function [Ansx, Ansy] = solve(ng, BB, delta)    
        %Guess
        y0 = Inf;
        y1 = [fsigmacrudexx(ng, BB); sigmaxy(ng, BB)];

        %Solve for sigmaemt
        while sum((abs(y1 - y0) ./ y1) >= delta)
            y0 = y1;
            y1(1) = sqrt(fiter1(ng, BB, y1)); %WATCH for sqrt
            y1(2) = fiter2(ng, BB, y1);
        end    

        Ansx = y1(1);
        Ansy = y1(2);
    end
%% Results
onesmatrix = ones(size(eta .* B));
[sigmaxxemt, sigmaxyemt] = arrayfun(@(nn, BB) solve(nn, BB, 1e-6), onesmatrix .* eta, onesmatrix .* B);

end


%% Helper Function
function Ans = sigmasimple(n, mu, nshift)

if nargin < 3 || isempty(nshift), nshift = 0; end

Ans = abs(n - nshift).*mu;
end
