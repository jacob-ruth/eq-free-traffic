function graphLiftN(hways, v0, evecs, evals, eps, oldData, t)
options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
h = 1.2;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 60;           % number of cars

posns = cumsum(hways(1:60));
system = [posns ; optimalVelocity(hways, v0)];
[tevo,evo] = ode45(@microsystem,[0 t],system, options,v0);
tsteps = linspace(0, t, t);
evoTPositions = interp1(tevo, evo, tsteps/2);


figure;
leslieAnn = getKnopePerkins(evoTPositions);
scatter(tsteps, leslieAnn(1,:), '.b');
figure;
scatter(tsteps, leslieAnn(2,:), '.r');


    function knopePerkins = getKnopePerkins(evo)
        knopePerkins = zeros(2,size(evo,1));
        for i = 1:size(evo,1)
            knopePerkins(:, i) = diffMapRestrict(getHeadways(evo(i,1:numCars)'), evals, evecs, oldData, eps);
        end
        
    end

% plot([0 t], [target target], 'g-');
%helpers that should be in a different file
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

%% function to circshift max to beginning
    function c = shiftMax(hways, center)
        if(nargin>1)
            shift = center + 1;
        else
            shift = 1;
        end
        [~, maxH] = max(hways,[],1);  % locate the max headway for each data point
        c = zeros(size(hways));
        
        % align all of the headways with the max in the front
        for iCar = 1:length(maxH)
            c(:,iCar) = circshift(hways(:,iCar), [-maxH(iCar)+shift,0]);
        end
    end
    function c = shiftMaxVel(hways, center)
        if(nargin>1)
            shift = center + 1;
        else
            shift = 1;
        end
        [~, maxH] = max(hways(1:numCars,:),[],1);  % locate the max headway for each data point
        c = zeros(size(hways));
        
        % align all of the headways with the max in the front
        for iCar = 1:length(maxH)
            c(1:numCars,iCar) = circshift(hways(1:numCars,iCar), [-maxH(iCar)+shift,0]);
            c(numCars+1:end, iCar) = circshift(hways(numCars+1:end,iCar), [-maxH(iCar)+shift,0]);
        end
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