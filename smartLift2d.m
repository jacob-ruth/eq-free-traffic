function lifted = smartLift2d(newVal, evec, eval, eps,v0, oldData, center)

options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
h = 1.2;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 60;           % number of cars
newVal = reshape(newVal, 1, length(newVal));        % newVal is a row vector
evoTime = 10;

tolerance = 1e-9;

distmult = 3;

foptions = optimset('TolFun',1e-10, 'ObjectiveLimit', 1e-6, 'Display', 'iter','Jacobian','on'); % fsolve options

difs = evec - repmat(newVal,size(evec,1),1);
allDists = sqrt(diag(difs * difs'));
[closestdist, closeidx] = min(allDists);
[closestprof,closestvals] = evolveRestrict([cumsum(oldData(:,closeidx)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx)),len),v0)]);

displacement = newVal - closestvals;    %vector from target to close (hopefully) point

pts = zeros(2,2);
perp = [-displacement(2) displacement(1)];
pts(1,:) = newVal + distmult*displacement - perp;
pts(2,:) = newVal + distmult*displacement + perp;

difs1 = evec - repmat(pts(1,:),size(evec,1),1);
allDists1 = sqrt(diag(difs1 * difs1'));
[closestdist1, closeidx1] = min(allDists1);
[closestprof1,closestvals1] = evolveRestrict([cumsum(oldData(:,closeidx1)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx1)),len),v0)]);

difs2 = evec - repmat(pts(2,:),size(evec,1),1);
allDists2 = sqrt(diag(difs2 * difs2'));
[closestdist2, closeidx2] = min(allDists2);
[closestprof2,closestvals2] = evolveRestrict([cumsum(oldData(:,closeidx2)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx2)),len),v0)]);

triangle = zeros(2*numCars,3);                   % points (profiles) of triangle in columns
triangle(:,1) = closestprof;
triangle(:,2) = closestprof1;
triangle(:,3) = closestprof2;

vals = [closestvals ; closestvals1 ; closestvals2];     % Leslie-Ann pairs (in rows) surrounding newVal

first = true;
done = false;
iter = 1;
tic;
figure;
while((first || ~done) && iter <= 200)
    iter = iter + 1;
    dists = [closestdist closestdist1 closestdist2];
% %     weights = (1/(2*sum(dists)))*(ones(length(dists),1)*sum(dists) - dists');
%     weights = ones(length(dists),1)/length(dists);
%     newprof = triangle*weights;
%     [dropprof, newpoint] = evolveRestrict(newprof);
    % find longest side length
    sides = squareform(pdist(vals));
    [~, location] = max(sides(:));
    [i,j] = ind2sub(size(sides),location);
    % split longest side
    newprof = (triangle(:,i) + triangle(:,j))/2;
    [dropprof, newpoint] = evolveRestrict(newprof);
   
    clf;
    hold on;
    scatter(vals(:,1),vals(:,2),'b*');
    scatter(newpoint(1),newpoint(2),'c*');
    scatter(newVal(1),newVal(2),'r*');
    hold off;
    pause;
    
    if(norm(newpoint - newVal) < tolerance)
        done = true;
        lifted = dropprof;
    elseif(~first && norm(oldVal - newpoint) < 1e-15)
        done = true;
        lifted = lineartriangle(vals,triangle);
    else
        if(PointInTriangle(newVal,[vals([1 2],:) ; newpoint]))
            triangle = [triangle(:,1:2), newprof];
            vals = [vals(1:2,:); newpoint];
        elseif(PointInTriangle(newVal,[vals([1 3],:) ; newpoint]))
            triangle = [triangle(:,[1 3]), newprof];
            vals = [vals([1 3],:); newpoint];
        elseif(PointInTriangle(newVal,[vals([2 3],:) ; newpoint]))
            triangle = [triangle(:,[2 3]), newprof];
            vals = [vals([2 3],:); newpoint];
        else
            warning('NOT IN A TRIANGLE (WELL TECHNICALLY ITS IN LOTS OF TRIANGLES BUT NONE OF OURS)');
        end
        
        oldVal = newpoint;
        first = false;
    end
end
disp(iter);
toc;

    function l = lineartriangle(vals,triangle)
        iterCount = 0;
        [lowpt,lowerPreDrop,highpt,higherPreDrop] = getBounds(vals,triangle);
        
        while(iterCount == 0 || (norm(val - newVal) > tolerance && iterCount < 50))
            lWeight = 0.5;                                                          % weight given to lower profile
            guessPrevo = lowerPreDrop * lWeight + higherPreDrop * (1 - lWeight);    % linear combination of upper and lower
            [guessEvolved,val] = evolveRestrict(guessPrevo);
            iterCount = iterCount + 1
            
%             clf;
%             hold on;
%             scatter(newVal(1),newVal(2),300,'r.');
%             scatter([lowpt(1),highpt(1)],[lowpt(2),highpt(2)],300,'b.');
%             scatter(val(1),val(2),'k*');
%             hold off;
%             pause;
            
            % reset the uppper or lower limit
            [lowpt,lowerPreDrop,highpt,higherPreDrop] =...
                getBounds([lowpt ; highpt ; val],[lowerPreDrop higherPreDrop guessPrevo]);
        end
        l = guessEvolved;
        fprintf('\t off by %d \n', norm(newVal-val));
    end

    function [pt1, prof1, pt2, prof2] = getBounds(pts,profs)
        distvecs = repmat(newVal,3,1) - pts;
        dist = sqrt(diag(distvecs * distvecs'));
        [dist,didx] = sort(dist);
        pts = pts(didx,:);
        profs = profs(:,didx);
        pt1 = pts(1,:);
        prof1 = profs(:,1);
        
        if((newVal-pt1)*(pts(2,:)-pt1)' >= 0)     % angle is acute ==> pts(2,:) is a valid choice
            pt2 = pts(2,:);
            prof2 = profs(:,2);
        else                                        % angle is obtuse ==> choose other point
            pt2 = pts(3,:);
            prof2 = profs(:,3);
        end
    end

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
        
        eventFunction = @(t,prof)loopEvent(t,prof,startAngle);
        
        options.events = eventFunction;
        
        [~,evolved] = ode45(@microsystem,[0 .5],cars,options,[v0 len h]);
%         while(abs(currentAngle - startAngle) > 0.05)                    % evolve until we're back at the starting angle
%             [~,evolved] = ode45(@microsystem,[0 .5],cars,options,[v0 len h]);
%             hways = getHeadways(evolved(end,1:numCars)', len);
%             restricted = diffMapRestrict(hways, eval, evec, oldData, eps);
%             xD = restricted(1) - center(1);
%             yD = restricted(2) - center(2);
%             currentAngle = atan2(yD,xD);
%             cars = evolved(end, :)';
%         end
        cars = evolved(end,:)';
        restricted = diffMapRestrict(getHeadways(evolved(end,1:numCars)',len),...
            eval, evec, oldData, eps);
        l = cars;
        val = restricted';
    end

    function [dif,isTerminal,direction] = loopEvent(~,prof,startAngle)
        isTerminal = 1;
        direction = 0;
        hways = getHeadways(prof(1:numCars)', len);
        restricted = diffMapRestrict(hways, eval, evec, oldData, eps);
        xD = restricted(1) - center(1);
        yD = restricted(2) - center(2);
        currentAngle = atan2(yD,xD);
        dif = currentAngle - startAngle;
    end

    function res = SameSide(p1,p2, a,b)
        a = [a 0];
        b = [b 0];
        p1 = [p1 0];
        p2 = [p2 0];
        cp1 = cross(b-a, p1-a);
        cp2 = cross(b-a, p2-a);
        if(dot(cp1, cp2) >= 0)
            res = true;
        else
            res = false;
        end
    end

    function res = PointInTriangle(p, t)
        a = t(1,:);
        b = t(2,:);
        c = t(3,:);
        if (SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b))
            res = true;
        else
            res = false;
        end
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