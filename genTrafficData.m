%% Generates data for traffic
% dataPoints-  number of different simulations to run
% dataPointsPer- the number of different sample points, evenly spaced, for
%                   each simulation
% v0-           the v0 to run each simulation
% muMin-        min perturbation for sin initial condition
% muMax-        max perturnbation for sin initial condition
%               actual mu will be evenly distributed between muMin and
%               muMax
% tMax-         the final time to stop the simulation at
% dataFile-     output file name, should be .mat
% plot-         bool, if it should be plot or not
function genTrafficData(dataPoints,dataPointsPer, v0, muMin, muMax,tMax, dataFile, plot)
h = 1.2;
len = 60;
numCars = 60;
options = odeset('AbsTol',10^-8,'RelTol',10^-8);
trafficOutput2 = zeros(2*numCars, dataPoints*dataPointsPer);
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
    scatter(t, std(getHeadways(allTime_1(:,1:numCars)')), '.')    
    figure;
    mu = muMax;
    finalTime = tMax;
    for i = 1:(numCars)
        cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars(i+numCars) = optimalVelocity(len/numCars, v0);
    end
    [t,allTime_1] = ode45(@microsystem,[0 finalTime],cars, options, v0);
    scatter(t, std(getHeadways(allTime_1(:,1:numCars)')), '.')
    title('max \mu')
    return;
end

for p = 1:dataPoints
    mu = muMin + (p/dataPoints)*(muMax - muMin);
    finalTime = tMax;
    for i = 1:(numCars)
        cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars(i+numCars) = optimalVelocity(len/numCars, v0);
    end
    [t,allTime_1] = ode45(@microsystem,[0 finalTime],cars, options, v0);
    tPoints = linspace(0, finalTime,dataPointsPer);
    points = interp1(t,allTime_1,tPoints)';
    trafficOutput2(:,((p - 1)*dataPointsPer + 1):((p)*dataPointsPer)) = points;
end;
save(dataFile,'trafficOutput2', 'v0', 'muMin', 'muMax', 'tMin', 'tMax');



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