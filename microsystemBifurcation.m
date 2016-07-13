function microsystemBifurcation()
h = 1.2;
len = 60;
numCars = 60;
mu = .1;

%% initialize car positions and velocities
cars_1 = zeros(2*numCars, 1);
cars_2 = zeros(2*numCars, 1);

origv01 = 0.91;
origv02 = 0.9;
v0_base1 = origv01;
v0_base2 = origv02;

for i = 1:numCars
    cars_1(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
    cars_1(i+numCars) = optimalVelocity(len/numCars, v0_base1);
    
    cars_2(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
    cars_2(i+numCars) = optimalVelocity(len/numCars, v0_base2);
end

%% save new reference state (if you want)
%{
tic;
[t1,allTime_1] = ode45(@microsystem,[0 finalTime],cars_1, options, v0_base1);
ref_1 = allTime_1(end,:)';
[t2,allTime_2] = ode45(@microsystem,[0 finalTime],cars_2, options, v0_base2);
ref_2 = allTime_2(end,:)';
toc;
save('refStats.mat','ref_1','ref_2');
%}

%% or you can load them instead
load('refStats919.mat','ref_1','ref_2');
hways1 = getHeadways(ref_1(1:numCars));
hways2 = getHeadways(ref_2(1:numCars));

%% Bifurcation on micro
steps = 110;
stepSize = 0.025;

cguess = -0.8581;

foptions = optimset( 'TolFun',1e-15,'TolX',1e-15, 'Display', 'off');

[~, max1] = max(hways1);
[~, max2] = max(hways2);
hways1 = circshift(hways1, [-max1  + 30, 0]);    % shift so max headway is in the middle-ish
hways2 =circshift (hways2, [-max2 + 30, 0]);
discreteScaling = 1;

% add more discretization points for smoother curve
hways1 = interp1(1:60,hways1,...
    linspace(1,numCars,numCars * discreteScaling))';
hways2 = interp1(1:60,hways2,...
    linspace(1,numCars,numCars * discreteScaling))';

sys1 = [hways1; cguess; 0; v0_base1];
sys2 = [hways2; cguess; 0; v0_base2];

bif  = zeros(numCars + 3, steps);

for iMic = 1:steps
    w = sys2 - sys1;        % tangent vector
    change =  stepSize*(w/norm(w));
    newGuess = sys2 + change;       % initial guess for newton solver
    lastGuess = sys2;               % reference state for phase condition
    [u , ~] = fsolve(@(sys2)FW(sys2, lastGuess, numCars, len, 1/1.7, w, newGuess),...
        newGuess, foptions);
    r = u(1:discreteScaling: end  - 3);
    bif(:,iMic) = [r; u(end - 2: end)];
    
    % set variables for next iteration
    sys1 = sys2;
    sys2 = u;
end

%% function to minimize with fsolve
% var       - state to vary in order to minimize fw
% ref       - reference state for phase condition
% nCars     - number of cars
% len       - track length
% tau       - momentum constant
% W         - tangent vector
% initGuess - first guess (point on tangent vector)
    function fw = FW(var, ref, nCars, len, tau, W, initGuess)
        v0 = var(end);
        fw = zeros(nCars*discreteScaling+3,1);
        fw(1:nCars*discreteScaling+2) = microFJ2(var(1:end-1), ref, nCars,...
            discreteScaling*nCars , len, v0, tau);
        fw(end) = W' * (var-initGuess);
    end

%% optimal velocity function given in paper
%    h         -  a parameter to represent the optimal velocity of the
%               	car
%    headway   - the distance between this car and the car ahead of it
%    v0        - ideal goal speed of each driver
%    returns:
%    v         - the optimal velocity of this car, who will speed or slow to try
%                   to meet it
    function v = optimalVelocity(headway,v0)
        v = v0 * (tanh(headway - h) + tanh(h));
    end

%% function to calculate headways of car vector
%  v   - column vector of the cars' positions
%  returns:
%  hways - column vector of cars' headways
    function hways = getHeadways(v)
        futureCars = circshift(v,[-1,0]);
        hways = mod(futureCars - v, len);
    end
end