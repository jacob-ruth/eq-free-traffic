function l = diffMapLift(newVal, evec, eval, eps, oldData)

fprintf('lifting to %d \n',newVal);

[~,closeidx] = min(abs(evec-newVal));
close = oldData(:,closeidx);

soptions = saoptimset('ObjectiveLimit',1e-6, 'AnnealingFcn', @changeWave, 'InitialTemperature', 0.75);
soptions.ReannealInterval = 100;

[l,fval] = simulannealbnd(@(x)toMin(x,newVal,evec,eval,eps,oldData),close,...
    zeros(length(close),1),60*ones(length(close),1),...
    soptions);
    %minimization function%
    function f = toMin(x,target,evec,eval,eps,old)
        lambda = 1;
        f = lambda * norm(diffMapRestrict(x,eval,evec,old,eps) - target);
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

end