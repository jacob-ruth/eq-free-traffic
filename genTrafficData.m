function genTrafficData(dataPoints, v0, muMin, muMax,tMin, tMax, dataFile, plot)
h = 1.2;
len = 60;
numCars = 60;
options = odeset('AbsTol',10^-8,'RelTol',10^-8);
buildSteps = dataPoints;
trafficOutput = zeros(2*numCars, buildSteps);
cars = zeros(numCars*2, 1);
if(plot)
    figure;
    mu = muMin;
    finalTime = tMax;
    for i = 1:(numCars)
        cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars(i+numCars) = optimalVelocity(len/numCars, v0);
    end
    [t,allTime_1] = ode45(@microsystem,[0 finalTime],cars, options, v0);
    %trafficOutput(:, p) = allTime_1(end,: )';
    scatter(t, std(getHeadways(allTime_1(:,1:numCars)')), '.')
    title('min \mu')

    figure;
    mu = muMax;
    finalTime = tMax;
    for i = 1:(numCars)
        cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars(i+numCars) = optimalVelocity(len/numCars, v0);
    end
    [t,allTime_1] = ode45(@microsystem,[0 finalTime],cars, options, v0);
    %trafficOutput(:, p) = allTime_1(end,: )';
    scatter(t, std(getHeadways(allTime_1(:,1:numCars)')), '.')
    title('max \mu')
end

for p = 1: buildSteps
    mu = muMin + rand*(muMax - muMin);
    finalTime = tMin + rand*(tMax - tMin);
    for i = 1:(numCars)
        cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars(i+numCars) = optimalVelocity(len/numCars, v0);
    end
    [t,allTime_1] = ode45(@microsystem,[0 finalTime],cars, options, v0);
    trafficOutput(:, p) = allTime_1(end,: )';
end;
save(dataFile,'trafficOutput', 'v0', 'muMin', 'muMax', 'tMin', 'tMax');



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

%% ODE that governs individual cars
% Governs the movement of the individual cars (microvariables)
% ~         - dummy parameter for time, to allow use in ode45
% params    - a column vector of the distribution of the cars and their
%               velocities, of size 2*numCars.  The position of car i
%               and its velocity are given at num params(i),
%               params(i + numcars)
%
    function u = microsystem(~,colCars, v0)
        invT = 1.7;
        headways = getHeadways(colCars(1:numCars));
        
        u = zeros(2*numCars,1);
        u(1:numCars,1) = colCars(numCars+1:2*numCars,1);
        u(numCars+1:2*numCars,1) = invT*(optimalVelocity(headways,v0) - colCars(numCars+1:2*numCars,1));
    end
end