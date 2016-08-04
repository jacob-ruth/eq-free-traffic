%% smartLift lifts from a diffusion map with 1 coordinate embedding
% newVal - target coordinate
% evec - diffusion map eigenvector
% eval - diffusion map eigenvalue
% eps - diffusion map epsilon
% v0 - optimal velocity parameter
% oldData - data used to build the diffusion map
% alignTo - location to align the max headways
% h - optimal velocity parameter
% len - length of the ring road
% Returns:
% l - the lifted profile
function  l = smartLift(newVal, evec, eval, eps,v0, oldData, alignTo, h, len)
options = odeset('AbsTol',10^-8,'RelTol',10^-8);    % ODE 45 options
numCars = size(oldData,1);                          % number of cars

preDrop = newVal;                                   
[~, valDropped] = dropDown(preDrop);                % closest to newVal evolved
dropDiff = valDropped - newVal;                     % distance from target and evolved closest
mult = 2;                                           % scaling factor for the linear combination
        
if( valDropped > newVal)                            % initalize the lower and upper limits
    higherP = preDrop;
    lowerP  = preDrop - mult*dropDiff;
else
    lowerP = preDrop;
    higherP  = preDrop - mult*dropDiff;
end

[~, valDroppedLower, lowerPreDrop] = dropDown(lowerP);      % find the evolution of the lower and upper limits
[~, valDroppedHigher, higherPreDrop] = dropDown(higherP);

if(valDroppedHigher < newVal || valDroppedLower > newVal)   % Check that newVal is between lower and higher
    warning('NewVal out of range! Increase mult');
end

% find the best lifted profile iteratively
val = 0;
iterCount = 0;
while(iterCount == 0 || (norm(val - newVal) > 1e-9 && iterCount < 20))
    lWeight = 0.5;                                                          % weight given to lower profile
    guessState = lowerPreDrop * lWeight + higherPreDrop * (1 - lWeight);    % linear combination of upper and lower
    guessPosns = cumsum(guessState);                                        % find positions 
    guessPrevo = [guessPosns ; optimalVelocity(h, guessState, v0)];         % inital state for ODE45
    [gtevo,gevo] = ode45(@microsystem,[0 10],guessPrevo, options,[v0 len h]);   % evolve the guess
    guessEvolved = getHeadways(gevo(end, 1:numCars)', len);                     
    guessShifted = alignMax(guessEvolved, alignTo, false);                  % align the max headway  
    val = diffMapRestrict(guessShifted(1:numCars),eval,evec,oldData,eps);   % restrict the evolved guess
    iterCount = iterCount + 1;
    
    % reset the uppper or lower limit
    if(val < newVal)
        lowerPreDrop = guessState;
    else 
        higherPreDrop = guessState;
    end
    
end

l = gevo(end, :);                 % return the evolved profile

    %% dropDown finds the closest profile to preDropVal, evolves and 
    % restricts it, and returns the resulting profile and embedding
    % preDropVal - target embedding coordinate
    % Returns:
    % l - the evolved profile with headways aligned
    % val - the restriction of l
    % preDrop - the closest value to preDropVal
    function [l, val, preDrop] = dropDown(preDropVal)
        [~, id] = min(abs(evec - preDropVal));          % find the closest value to preDropVal
        neighbor = oldData(:, id);                      % the closest profile
        if(nargout == 3)
            preDrop = neighbor;     
        end
        liftedPosns = cumsum(neighbor);                 % find the positions
        liftedGuess = [liftedPosns ; optimalVelocity(h, neighbor, v0)];     % initial state for ODE45
        [~,evo] = ode45(@microsystem,[0 10],liftedGuess, options,[v0 len h]);   % evolve the profile
        evolvedHways = getHeadways(evo(end, 1:numCars)', len);  
        shifted = alignMax(evolvedHways, alignTo, false);             	 % align the max headway
        val = diffMapRestrict(shifted,eval,evec,oldData,eps);           % restrict the evolved profile
        l = shifted;
    end


end