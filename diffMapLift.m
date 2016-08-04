function l = diffMapLift(newVal, evec, eval, eps, oldData)
numCars = 60;
[diffs, idx] = sort(abs(evec - newVal));
numVals = 4;
neighbors = oldData(:, idx(1:numVals));
realDifs = diffs(1:numVals);
neighborDifs = exp(-(realDifs)/(1e-5));
f =  @(x, problem)linearComb(x, problem,neighbors, neighborDifs);
soptions = saoptimset('ObjectiveLimit',5e-5, 'AnnealingFcn', f,...
    'InitialTemperature', 0.5, 'Display', 'off','MaxIter',1000);
soptions.ReannealInterval = 100;
close = neighbors(:,1);
[l,fval] = simulannealbnd(@(x)toMin(x,newVal,evec,eval,eps,oldData),close,...
    zeros(numCars,1),60*ones(numCars,1),...
    soptions);  
% l = close;
% figure;
% hold on;
% scatter(1:60, closest(1:60), 'rx');
% scatter(1:60, l(1:60), 'bo');
% scatter(1:60, close(1:60), 'g.');
% % 

% fprintf('Lifting finished after %i iterations, with fval of %d, and difference from close of %d\n', output.iterations, fval, norm(l - close));

%     %minimization function%
    function f = toMin(x,target,evec,eval,eps,old)
        lambda = 1;
        f = lambda * norm(diffMapRestrict(x(1:numCars),eval,evec,old,eps) - target);
    end
    %changes wave based on distance from average
    function newx = changeWave(optimValues, problem)
        currPositions = optimValues.x;
        avgHway = sum(currPositions) / length(currPositions);
        hwayDif = currPositions - avgHway;
        a = rand() - 0.5;
        b = (a * optimValues.temperature);
        newx = currPositions + hwayDif.*b;
    end

    %changes wave based on curvature
    function newx = heatDiffuse(optimValues, problem)
        old = optimValues.x;
        LD2 = fourdif(length(old),2)*2*pi/length(old);
        D2 = LD2 * old;
        r = -0.05 + rand()/10;
        new = old + r*D2;
        new = smooth(new, 'sgolay');
        newx = new - ((sum(new) - sum(old))/ sum(old));
    end

    function newx = linearComb(optimValues, problem, neighbors, neighDiffs)
       hways = optimValues.x(1:numCars);
       vels = optimValues.x(numCars+1: end);
       weights = rand(length(neighDiffs), 1).*neighDiffs;
       weights = weights/sum(weights);
       newA = neighbors*(weights);
       newx = hways.*(1 - optimValues.temperature(1:60)) + newA.*optimValues.temperature(1:60);     
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
    
end