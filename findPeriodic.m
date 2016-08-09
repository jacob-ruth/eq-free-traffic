%% finds the evolved state that has gone around the ring road exactly once
% cars  - 2*numCars vector of positions (NOT HEADWAYS) and velocities
% eval  - eigenvalues from diffusion map
% evec  - eigenvectors from diffusion map in columns
% eps   - epsilon from diffusion map
% v0    - v0 (optimal velocity) parameter
%
% returns: profile of positions and velocities arising from evolving cars
%          for one period
function [looped, tangent] = findPeriodic(cars, eval, evec, oldData, eps, v0, tang, started)
len = 60;
numCars = length(cars)/2;
h = 1.2;

if(nargin < 7)
    hways1 = getHeadways(cars(1:numCars), len);
    startRestrict = diffMapRestrict(hways1, eval, evec, oldData, eps);     % find the starting coordinate
    tangentLine = finiteDifference(cars);
    minLoopTime = 65;
else
    tangentLine = tang;
    startRestrict = started;
    minLoopTime = 0;
end

eventFunction = @(t,prof)loopEvent(t,prof, [tangentLine startRestrict]);

opts = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
opts = odeset(opts,'Events',eventFunction);

[~,evolve,te] = ode45(@(t,y)microsystem(t,y,[v0 len h]), [0 200], cars, opts);

if(length(te) < 1)
    fprintf('No periodic solution found');
end

looped = evolve(end, :)';
if(nargout > 1)
    tangent = tangentLine;
end

    function diff = finiteDifference(cars)
        options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
        hwaysStart = getHeadways(cars(1:numCars), len);
        start = diffMapRestrict(hwaysStart, eval, evec, oldData, eps);
        [~,evol] = ode45(@microsystem,[0 .01],cars,options,[v0 len h]);
        hwaysEnd = getHeadways(evol(end,1:numCars)', len);
        final = diffMapRestrict(hwaysEnd, eval, evec, oldData, eps);
        diff = start - final;
    end

    function [dif,isTerminal,direction] = loopEvent(t,prof, param)
        tan = param(:,1);
        start = param(:,2);
        
        if(t > minLoopTime)
            isTerminal = 1;
            hways = getHeadways(prof(1:numCars), len);
            restricted = diffMapRestrict(hways, eval, evec, oldData, eps);
            current = restricted - start;
            dif = dot(current, tan);
        else
            isTerminal = 0;
            dif = .1;
        end
        direction = -1;
    end
end