function l = diffMapLiftAnneal(newVal, evec, eval, eps, oldData)

%% initialize simulated annealing
closeidxs = dsearchn(evec,newVal');     % find the initial guess
closeidx = closeidxs(1);
x = oldData(:,closeidx);

optimal = x;                                            % keep track of a global optimal solution
oldCost = toMin(x, newVal, evec, eval, eps, oldData);	% initialize the cost
optimalCost = oldCost;                                  % initialize the optimal cost
temp = 1.0;                                             % starting temperature
tempMin = 0.00001;                                      % ending temeprature
alpha = 0.9;
tolerance = 1.5*10^(-7);

soptions = saoptimset('ObjectiveLimit',1e-5,'AnnealingFcn',@annealingboltz,'InitialTemperature',1e-3);
soptions.ReannealInterval = 50;

%% run simulated annealing
%{
while(temp > tempMin && optimalCost > tolerance)
    i = 1;
    while(i <= 100)
        newX = getNeighbor(x);             % find the new neighbor
        %plot(1:1:60, newX); drawnow; pause;
        newCost = toMin(newX, newVal, evec, eval, eps, oldData);	% calculate the cost at the new location
        ap = acceptanceProbability(oldCost, newCost, temp);         % find the acceptance probability
        if (ap > rand())                                         	% move to this neighbor based on the acceptance probability
            if(newCost < optimalCost)                               % if this is a new best solution, save it as the current optimal
                optimal = newX;
                optimalCost = newCost;
                fprintf('\t Found a new optimal cost: %e \n', optimalCost);
            end
            x = newX;
            oldCost = newCost;
        end
        i = i + 1;
    end
    temp = temp*alpha;                      % decrease the temperature
end

l = optimal;                                % return the best solution found
%}

% function to optimize
    function f = toMin(x,target,evec,eval,eps,old)
        lambda = 1;
        f = lambda * norm(diffMapRestrict(x,eval,evec,old,eps)-target);
        f = f + abs(sum(x) - 60);
    end

% probability to accept a new solution
    function ap = acceptanceProbability(old, new, t)
        ap = exp((old-new)/t);
    end

% find the neighbor for simulated annealing
    function new = getNeighbor(old)
        LD2 = fourdif(length(old),2)*2*pi/length(old);
        D2 = LD2 * old;
        r = -0.05 + rand()/10;
        new = old + r*D2;
        new = smooth(new, 'sgolay');
        new = new - ((sum(new) - sum(old))/ sum(old));
    end
end