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
function [lifted, idxMin] = newLift(newVal, evec, eval, eps,v0, oldData)
lsqOptions = optimset('Display','off');
lsqOptions.FunctionTolerance = 1e-8;
lsqOptions.StepTolerance = 1e-9;
%% constants
h = 2.4;                % optimal velocity parameter
len = 60;               % length of the ring road

minPoints = 10; % take linear combination of at least 10 profiles
%radius = 0.00004;
%radiusIncrement = 0.00001;

%% find closest datapoints

validPoints = 10;

newDist = pdist2(newVal',evec)';

[~, index] = sort(newDist,'ascend');
idxMin = index(1:minPoints);
closestPoints = oldData(:, idxMin); % use closest profiles for linear combination
validPoints = minPoints;
idxMin = sort(idxMin);

%{
newKernel = exp(-newDist.^2/eps^2); % calculate weights from previous data and new point
while validPoints < minPoints
    radius = radius + radiusIncrement;
    minKernel = exp(-1 * radius ^ 2 / eps); % why only dividing by epsilon, not epsilon squared?
    withinRange = (newKernel > minKernel) & (newKernel ~= 1); % count points that have at least this weight, why newKernel ~= 1? If it's already in the set, that should be good
    validPoints = sum(withinRange);
end

% closestEvecs = evec(withinRange, :);
closestPoints = oldData(:, withinRange); % use closest profiles for linear combination
%}
% solve for best coefficients for linear combination
liftedCoeffs = lsqnonlin(@liftingEquations, 1/validPoints*ones(validPoints, 1), zeros(validPoints, 1), ones(validPoints,1), lsqOptions);
%liftedCoeffs = fsolve(@liftingEquations, zeros(validPoints, 1));
liftedHeadway = closestPoints * liftedCoeffs; % create new lifted profile
lifted = [hwayToPos(liftedHeadway) ; optimalVelocity(h, liftedHeadway, v0)];

%disp('error');
%disp(norm((newVal-diffMapRestrictAlt(liftedHeadway,eval,evec,oldData,eps,1))./newVal)*100);

%% sets up and solves a system of equations

    function out = liftingEquations(coeffs)
        liftGuess = closestPoints * coeffs;
        restrictGuess = diffMapRestrictAlt(liftGuess,eval,evec,oldData,eps,1);
        %out = [restrictGuess - newVal; sum(coeffs) - 1; abs(coeffs) - coeffs];
        out = [norm(restrictGuess) - norm(newVal);  atan2(restrictGuess(2),restrictGuess(1)) - ...
            atan2(newVal(2), newVal(1)) ; sum(coeffs) - 1; zeros(validPoints - 2, 1)]; 
    end

end