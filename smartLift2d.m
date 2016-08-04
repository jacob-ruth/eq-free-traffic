function [closestdist, l] = smartLift2d(newVal, evec, eval, eps,v0, oldData, center)
options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
h = 1.2;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 60;           % number of cars
newVal = reshape(newVal, 1, length(newVal));        % newVal is a row vector
evoTime = 10;

distmult = 2;

foptions = optimset('TolFun',1e-10, 'ObjectiveLimit', 1e-6, 'Display', 'iter','Jacobian','on'); % fsolve options

difs = evec - repmat(newVal,size(evec,1),1);
allDists = sqrt(diag(difs * difs'));
[closestdist, closeidx] = min(allDists);
[closestprof,closestvals] = evolveRestrict([cumsum(oldData(:,closeidx)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx)),len),v0)]);

displacement = newVal - closestvals;    %vector from target to close (hopefully) point

triangle = zeros(2*numCars,3);                   % points of triangle in columns

triangle(:,1) = closestprof;

pts = zeros(2,2);
perp = [-displacement(2) displacement(1)];
pts(1,:) = newVal + distmult*displacement - perp;
pts(2,:) = newVal + distmult*displacement + perp;

figure;
hold on;
scatter([closestvals(1) ; pts(:,1)], [closestvals(2) ; pts(:,2)], 200, 'b.');
scatter(evec(closeidx,1),evec(closeidx,2),200,'g*');
scatter(newVal(1),newVal(2),200,'r.');
plot([newVal(1) newVal(1)+distmult*displacement(1)],[newVal(2) newVal(2) + distmult*displacement(2)],'g-');
plot([newVal(1)+distmult*displacement(1)+perp(1) newVal(1)+distmult*displacement(1)-perp(1)],[newVal(2)+distmult*displacement(2)+perp(2) newVal(2)+distmult*displacement(2)-perp(2)],'r-');

difs1 = evec - repmat(pts(1,:),size(evec,1),1);
allDists1 = sqrt(diag(difs1 * difs1'));
[closestdist1, closeidx1] = min(allDists1);
[closestprof1,closestvals1] = evolveRestrict([cumsum(oldData(:,closeidx1)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx1)),len),v0)]);

triangle(:,2) = closestprof1;

difs2 = evec - repmat(pts(2,:),size(evec,1),1);
allDists2 = sqrt(diag(difs2 * difs2'));
[closestdist2, closeidx2] = min(allDists2);
[closestprof2,closestvals2] = evolveRestrict([cumsum(oldData(:,closeidx2)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx2)),len),v0)]);

triangle(:,3) = closestprof2;

vals = [closestvals ; closestvals1 ; closestvals2];

scatter(evec([closeidx1 closeidx2],1),evec([closeidx1 closeidx2], 2), 200,'g*');

scatter(vals(:,1), vals(:,2), 200, 'k.');
scatter(newVal(1), newVal(2), 200, 'r.');
hold off;


first = true;
done = false;
tic;
while(first || ~done)
    first = false;
    numPts = size(evec,2) + 1;
    triangle = zeros(size(evec,2),numPts);       % points are in columns
    profiles = zeros(numCars, numPts);
    for iPt = 1:numPts
        [profiles(:,iPt), triangle(:,iPt)] = evolveRestrict(oldData(:,randi(size(oldData,2))));
    end
    
    if(inpolygon(newVal(1),newVal(2),triangle(1,:),triangle(2,:)))
        fprintf('WE DID A THING\n');
        done = true;
    else
        fprintf('we did not do a thing :(\n');
        done = false;
    end
end
toc;



    function [f,J] = toMin(val)
        lambda = 1;
        [~, dropped] = dropDown(val);
        f = lambda * norm(dropped - newVal);
        J = 1;
    end

    % cars - 2*numCars vector of positions (NOT HEADWAYS) and velocities
    function [l, val] = evolveRestrict(cars)
        hways = getHeadways(cars(1:numCars), len);
        start = diffMapRestrict(hways, eval, evec, oldData, eps);     % find the starting coordinate   
        xDist = start(1) - center(1);
        yDist = start(2) - center(2);
        startAngle = atan2(yDist,xDist);                                % find the angle to match
        currentAngle = 10000;
        while(abs(currentAngle - startAngle) > 0.05)                    % evolve until we're back at the starting angle
            [~,evolved] = ode45(@microsystem,[0 1],cars,options,[v0 len h]);
            hways = getHeadways(evolved(end,1:numCars)', len);
            restricted = diffMapRestrict(hways, eval, evec, oldData, eps);
            xD = restricted(1) - center(1);
            yD = restricted(2) - center(2);
            currentAngle = atan2(yD,xD);
            cars = evolved(end, :)';
        end
        l = cars;
        val = restricted';
    end

    function [l, val] = dropDown(preDropVal)
        [~, indx] = min();
        neighbor = oldData(:, indx);
        liftedPosns = cumsum(neighbor);
        liftedGuess = [liftedPosns ; neighborsV];
        [tevo,evo] = ode45(@microsystem,[0 evoTime],liftedGuess, options,v0);
        evolvedHways = getHeadways(evo(end, 1:numCars)', len);
        val = diffMapRestrict(evolvedHways(1:numCars),eval,evec,oldData,eps);
        l = evolvedHways;
    end

end