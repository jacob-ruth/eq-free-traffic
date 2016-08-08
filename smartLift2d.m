function lifted = smartLift2d(newVal, evec, eval, eps,v0, oldData)

h = 1.2;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 60;           % number of cars
newVal = reshape(newVal, 1, length(newVal));        % newVal is a row vector

tolerance = 5e-8;

distmult = 2;

difs = evec - repmat(newVal,size(evec,1),1);
allDists = sqrt(sum(difs.^2,2));
[~,sortidx] = sort(allDists);

radius = max(max(abs(evec(:,1)),max(abs(evec(:,2)))));
if(newVal(1)^2 + newVal(2)^2 >= radius^2)
    closeidx = sortidx(1);
    lifted = [cumsum(oldData(:,closeidx)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx)),len),v0)];
    fprintf('Lifting to value outside of data set \n');
    return;
end

vals = 10000*ones(2,3);
itri = 1;
tic;
while(itri == 1 || ~PointInTriangle(newVal,vals) )
    closeidx = sortidx(itri);
    [closestprof,closestvals] = evolveRestrict([cumsum(oldData(:,closeidx)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx)),len),v0)]);
    
    displacement = newVal - closestvals;    %vector from target to close (hopefully) point
    
    pts = zeros(2,2);
    perp = [-displacement(2) displacement(1)];
    pts(1,:) = newVal + distmult*displacement - perp;
    pts(2,:) = newVal + distmult*displacement + perp;
    
    difs1 = evec - repmat(pts(1,:),size(evec,1),1);
    allDists1 = sqrt(sum(difs1.^2,2));
    [~, closeidx1] = min(allDists1);
    [closestprof1,closestvals1] = evolveRestrict([cumsum(oldData(:,closeidx1)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx1)),len),v0)]);
    
    difs2 = evec - repmat(pts(2,:),size(evec,1),1);
    allDists2 = sqrt(sum(difs2.^2,2));
    [allDists2,~] = sort(allDists2);
    [~, closeidx2] = min(allDists2);
    [closestprof2,closestvals2] = evolveRestrict([cumsum(oldData(:,closeidx2)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx2)),len),v0)]);
    
    triangle = zeros(2*numCars,3);                   % points (profiles) of triangle in columns
    triangle(:,1) = closestprof;
    triangle(:,2) = closestprof1;
    triangle(:,3) = closestprof2;
    
    vals = [closestvals ; closestvals1 ; closestvals2];     % Leslie-Ann pairs (in rows) surrounding newVal

    itri = itri+1;
end
x = toc;
fprintf('Found a triangle in %d step(s) and %f seconds \n', itri-1, x);

first = true;
done = false;
notMedian = false;
iter = 1;
tic;
while((first || ~done) && iter <= 50)
    iter = iter + 1;
    
    % find longest side length
    sides = squareform(pdist(vals));
    [~, location] = max(sides(:));
    [i,j] = ind2sub(size(sides),location);
    if(notMedian)
        fprintf('Triangle stalled \n');
        i = mod(i+1,3) + 1;
        j = mod(j+1,3) + 1;
    end

    % split longest side
    weight = 0.5;
    newprof = weight*triangle(:,i) + (1-weight)*triangle(:,j);
    [dropprof, newpoint] = evolveRestrict(newprof);
    
    if(~first && norm(oldVal - newpoint) < 10^(-20))
        notMedian = true;
    else
        notMedian = false;
    end
%     clf;
%     hold on;
%     scatter(vals(:,1),vals(:,2),'b*');
%     scatter(newpoint(1),newpoint(2),'c*');
%     scatter(newVal(1),newVal(2),'r*');
%     if(~first)
%         scatter(oldVal(1),oldVal(2),'ko');
%     end
%     drawnow;
%     hold off;
    
    if(norm(newpoint - newVal) < tolerance)
        done = true;
        lifted = dropprof;
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
        
        first = false;
        oldVal = newpoint;
    end
end
x = toc;
fprintf('Triangle converged in %d step(s) and %f seconds \n', iter-1, x);

if(~done)
    fprintf('Triangle did not converge \n');
    dists = vals  - repmat(newVal, 3,1);
    distances = sqrt(sum(dists.^2,2));
    [~,sortidx] = sort(distances);
    lifted = triangle(:,sortidx(1));
end

% cars - 2*numCars vector of positions (NOT HEADWAYS) and velocities
    function [l, val] = evolveRestrict(cars)
        l = findPeriodic(cars,eval,evec,oldData,eps,v0);
        restricted = diffMapRestrict(getHeadways(l(1:numCars),len),...
            eval,evec,oldData,eps);
        val = restricted';
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

end