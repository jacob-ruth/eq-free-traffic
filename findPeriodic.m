function looped = findPeriodic(cars, eval, evec, oldData, eps, v0)
len = 60;
numCars = length(cars)/2;
h = 1.2;

hways1 = getHeadways(cars(1:numCars), len);
startRestrict = diffMapRestrict(hways1, eval, evec, oldData, eps);     % find the starting coordinate
tangentLine = finiteDifference(cars);

eventFunction = @(t,prof)loopEvent(t,prof, [tangentLine startRestrict]);

opts = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
opts = odeset(opts,'Events',eventFunction);

[~,evolve] = ode45(@(t,y)microsystem(t,y,[v0 len h]), [0 100], cars, opts);
looped = evolve(end, :)';


    function diff = finiteDifference(cars)
        options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
        hwaysStart = getHeadways(cars(1:numCars), len);
        start = diffMapRestrict(hwaysStart, eval, evec, oldData, eps);
        [~,evol] = ode45(@microsystem,[0 .01],cars,options,[v0 len h]);
        hwaysEnd = getHeadways(evol(end,1:numCars)', len);
        final = diffMapRestrict(hwaysEnd, eval, evec, oldData, eps);
        diff = start - final;
    end

    function [dif,isTerminal,direction] = loopEvent(~,prof, param)
        tan = param(:,1);
        start = param(:,2);
        
        isTerminal = 1;
        direction = -1;
        
        hways = getHeadways(prof(1:numCars), len);
        restricted = diffMapRestrict(hways, eval, evec, oldData, eps);
        current = restricted - start;
        dif = dot(current, tan);
    end
end