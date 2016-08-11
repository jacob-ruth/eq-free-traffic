%% Generates data for traffic
% dataPoints-  number of different simulations to run
% dataPointsPer- the number of different sample points, evenly spaced, for
%                   each simulation
% v0-           the v0 to run each simulation
% muMin-        min perturbation for sin initial condition
% muMax-        max perturnbation for sin initial condition
%               actual mu will be evenly distributed between muMin and
%               muMax
% tMin =        the minimum time to record the simulation at
% tMax-         the final time to stop the simulation at
% dataFile-     output file name, should be .mat
% plotIt-       bool, if it should be plot or not, if so there will only be
%                   1 simulation
function genTrafficData(dataPoints,dataPointsPer, v0, muMin, muMax,tMin, tMax, dataFile, plotIt)
h = 1.2;
len = 60;
numCars = 60;
options = odeset('AbsTol',10^-8,'RelTol',10^-8);
cars = zeros(numCars*2, 1);
timeStep = (tMax - tMin)/dataPointsPer;
trafficData = zeros(2*numCars, dataPoints*dataPointsPer);

% if plot, run one simulation
if(plotIt)
    mu = muMin + rand()*(muMax - muMin);
    for i = 1:(numCars)
        cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
        cars(i+numCars) = optimalVelocity(h, len/numCars, v0);
    end
    [t,allTime] = ode45(@microsystem,[0 tMin],cars, options, [v0 len h]);           % evolve to the minimum time
    sig = std(getHeadways(allTime(:, 1:numCars)', len));                % plot standard deviation
    figure;
    plot(t, sig');
    finalCars = allTime(end,:)';
    trafficData = finalCars;                                            % save the final simulation data
else
    % otherwise run dataPoints simulations and take dataPointsPer snapshots
    % of each simulation
    for iOuter=1:dataPoints
        trafficOutput2 = zeros(2*numCars, dataPoints);
        mu = muMin + rand()*(muMax - muMin);                % initialize the cars
        for i = 1:(numCars)
            cars(i) = (i-1) * len/numCars + mu*sin(2*pi*i/numCars);
            cars(i+numCars) = optimalVelocity(h, len/numCars, v0);
        end
        [~,allTime] = ode45(@microsystem,[0 tMin],cars, options, [v0 len h]);           % evolve to the minimum time
        finalCars = allTime(end,:)';
        trafficOutput2(:,1) = finalCars; 
        for p = 2:dataPointsPer                         % grab the snapshots
            [~,allTime] = ode45(@microsystem,[0 timeStep], finalCars, options, [v0 len h]); % evolve by the step
            finalCars = allTime(end,:)';
            trafficOutput2(:,p) = finalCars;
        end
        trafficData(:, (iOuter -1)*dataPointsPer + 1:iOuter*dataPointsPer) = trafficOutput2;
    end
end
vel = v0*ones(dataPoints*dataPointsPer, 1); % save the value of v0 used on this run

save(dataFile,'trafficData', 'vel');  % save the snapshots of the simulation 

end