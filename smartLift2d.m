%% two-dimensional lifting operator
% newVal    - 1x2 vector representing the desired Leslie and Ann values
% evec      - eigenvectors from diffusion map in columns
% eval      - eigenvalues from diffusion map
% eps       - epsilon from diffusion map
% v0        - v0 parameter value
% oldData   - data from which the diffusion map was constructed
%
% returns a 120 x 1 vector of car positions and velocities that restricts
% to close to newVal
function lifted = smartLift2d(newVal, evec, eval, eps,v0, oldData)

h = 1.2;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 60;           % number of cars
tolerance = 1e-12;       % lifting tolerance
distmult = 2;           % how far to go for the vertices of the triangle

newVal = reshape(newVal, 1, length(newVal));        % newVal is a row vector

%% calculate distance from each known data point to newVal
difs = evec - repmat(newVal,size(evec,1),1);
allDists = sqrt(sum(difs.^2,2));
[~,sortidx] = sort(allDists);

%% check if outside data set
radius = max(max(abs(evec(:,1)),max(abs(evec(:,2)))));
if(norm(newVal) >= radius)
    closeidx = sortidx(1);
    lifted = [cumsum(oldData(:,closeidx)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx)),len),v0)];
    fprintf('Lifting to value outside of data set \n');
    return;
end

%% loop through closest points until find one that makes a triangle around newVal
tic;
itri = 1;
vals = 1000*ones(3,2);
while(itri == 1 || ~PointInTriangle(newVal,vals) )
    [triangle, vals] = findTriangle(itri);
    itri = itri+1;
end
x = toc;
fprintf('\tFound a triangle in %d step(s) and %f seconds \n', itri-1, x);

%% split the longest side to converge on newVal
first = true;
done = false;
notMedian = false;
iter = 1;
tic;
while((first || ~done) && iter <= 500)
    iter = iter + 1;
    
    % if the triangle is stuck at the previous locations, create a new
    % initial triangle
    if(notMedian)
        fprintf('\tTriangle stalled \n');
        while(notMedian || ~PointInTriangle(newVal,vals) )
           notMedian = false;
           [triangle, vals] = findTriangle(itri);
           itri = itri+1;
        end
    end
    
    % find longest side length
    sides = squareform(pdist(vals));
    [~, location] = max(sides(:));
    [i,j] = ind2sub(size(sides),location);
    % split longest side
    weight = 0.5;
    newprof = weight*triangle(:,i) + (1-weight)*triangle(:,j);
    [dropprof, newpoint] = evolveRestrict(newprof);
    
    % check if the triangle is stuck at the same locations
    if(~first && norm(oldVal - newpoint) < 10^(-20))
        notMedian = true;
    else
        notMedian = false;
    end
    
    %{
    clf;
    hold on;
    scatter(vals(:,1),vals(:,2),'b*');
    scatter(newpoint(1),newpoint(2),'c*');
    scatter(newVal(1),newVal(2),'r*');
    if(~first)
        scatter(oldVal(1),oldVal(2),'ko');
    end
    drawnow;
    hold off;
    %} 
    
    % check if the triangle has converged yet
    if(norm(newpoint - newVal) < tolerance)
        done = true;
        lifted = dropprof;
    else % if not, reassign the vertices so newVal is still in the triangle
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
fprintf('\tTriangle converged in %d step(s) and %f seconds \n', iter-1, x);

% if the triangle doesn't converge, return the closest profile
if(~done)
    fprintf('\tTriangle did not converge \n');
    dists = vals  - repmat(newVal, 3,1);
    distances = sqrt(sum(dists.^2,2));
    [~,sortidx] = sort(distances);
    lifted = triangle(:,sortidx(1));
end

%% given a profile, evolve it for one period and restrict to Leslie/Ann
% cars - 2*numCars vector of positions (NOT HEADWAYS) and velocities
%
% returns: l    - the evolved profile
%          val  - the restricted Leslie/Ann values
    function [l, val] = evolveRestrict(cars)
        l = findPeriodic(cars,eval,evec,oldData,eps,v0);
        restricted = diffMapRestrict(getHeadways(l(1:numCars),len),...
            eval,evec,oldData,eps);
        val = restricted';
    end

%% detect whether p1 and p2 are on the same side of the line AB
% p1    - 2x1 vector representing the first point
% p2    - 2x1 vector representing the second point
% a     - 2x1 vector representing the one point on the line
% b     - 2x1 vector representing the second point defining the line
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


%% detects whether a point is in a given triangle
% p     - 2x1 vector representing the point to check
% t     - 3x2 matrix representing the triangle with the points in rows
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

%% function findTriangle tries to find a triangle that surrounds newVal
% itri - the index of the closest point to try for the first vertex 
% Returns:
% tri - the profiles of the triangle (positions and velocities)
% values - the diffusion map coordinates of the triangle
    function [tri, values] = findTriangle(itri)
        closeidx = sortidx(itri);   % find the first vertex at the itri closest position
        [closestprof,closestvals] = evolveRestrict([cumsum(oldData(:,closeidx)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx)),len),v0)]);
        
        displacement = newVal - closestvals;    %vector from target to close (hopefully) point
        
        pts = zeros(2,2);                           % find vertices that are (hopefully) on the other side of newVal
        perp = [-displacement(2) displacement(1)];
        pts(1,:) = newVal + distmult*displacement - distmult*perp;
        pts(2,:) = newVal + distmult*displacement + distmult*perp;
        
        difs1 = evec - repmat(pts(1,:),size(evec,1),1); % find the closest profile to the second vertex
        allDists1 = sqrt(sum(difs1.^2,2));
        [~, closeidx1] = min(allDists1);
        [closestprof1,closestvals1] = evolveRestrict([cumsum(oldData(:,closeidx1)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx1)),len),v0)]);
        
        difs2 = evec - repmat(pts(2,:),size(evec,1),1); % find the closest profile to the third vertex
        allDists2 = sqrt(sum(difs2.^2,2));
        [~, closeidx2] = min(allDists2);
        [closestprof2,closestvals2] = evolveRestrict([cumsum(oldData(:,closeidx2)) ; optimalVelocity(h,getHeadways(cumsum(oldData(:,closeidx2)),len),v0)]);
        
        % return the profiles and coordinates of the triangle
        tri = [closestprof closestprof1 closestprof2];            % profiles in columns
        values = [closestvals ; closestvals1 ; closestvals2];     % Leslie-Ann pairs (in rows) surrounding newVal

    end

end