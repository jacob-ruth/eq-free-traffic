function trafficBifurcation()
    h = 2.4;                % optimal velocity parameter
    len = 60;               % length of the ring road
    numCars = 30;           % number of cars

    tskip = 100;            % times for evolving
    delta = 1000;           
    stepSize = .0025;        % step size for the secant line approximation
    delSigma = 0.00001;     % delta sigma used for finite difference of F
    delv0 = 0.00001;        % delta v0 used for finite difference of F
    tolerance = 10^(-7);    % tolerance for Newton's method
    
    options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options

    % load the reference states if they have already been calculated
    load('microBif.mat', 'bif', 'vel');
    oldBif = bif;
    
    start = 105;
    ref_2 = oldBif(1:numCars, start+1);   % initial profiles
    ref_1 = oldBif(1:numCars, start);
    v0_base1 = vel( start);
    v0_base2 = vel(1 + start);
    
    %% initialize secant continuation
    steps = 100;                                % number of steps to take around the curve
    bif = zeros(2,steps);                       % array to hold the bifurcation values
    actuals = zeros(2,steps);

    figure; hold on; % draw the bifurcation diagram
    scatter(vel(start:end), std(oldBif(1:numCars,start:end)), 200, 'r.');
    
    ref_1 = hwayToPos(ref_1);                          % get reference state to position
    ref_2 = hwayToPos(ref_2);
    sigma_1 = std(getHeadways(ref_1(1:numCars)));      % initial sigma values for secant line approximation
    sigma_2 = std(getHeadways(ref_2(1:numCars)));
    
    %% pseudo arc length continuation
    for iEq=1:steps
        fprintf('Starting iteration %d of %d \n', iEq, steps);
        w = [sigma_2 - sigma_1 ; v0_base2 - v0_base1];          % slope of the secant line
        newGuess = [sigma_2; v0_base2] + stepSize *(w/norm(w)); % first guess on the secant line

        %% initialize Newton's method
        u = newGuess;
        first = true;                       % mimic a do-while loop                                
        k=1;                            	% Newton's method counter
        %% Newton and that other guy's method
        while(first ||(norm(invD*[f;neww])>tolerance && k < 20))
            first = false;
            fprintf('\tNewton iteration: %d \n', k);
            f = F(ref_2, u(1),u(2));                                % calculcate the function to zero
            neww = w(1)*(u(1)-newGuess(1)) + w(2)*(u(2) - newGuess(2));
            Df = jacobian(ref_2, u(1), u(2), w);                    % find the jacobian
            invD = Df^(-1);
            u = u - invD*[f;neww];                                  % perform the Newton step
            k = k + 1;
        end  
        
        %% alternate Newton's method using fsolve
        % u = fsolve(@(u)FW(u,ref_2,w,newGuess), newGuess); 

        bif(:,iEq) = u;                                             % save the new solution

        %% reset the values for the arc length continuation
        ref_1 = ref_2;
        sigma_1 = sigma_2;
        v0_base1 = v0_base2;
        v0_base2 = u(2);
        [sigma_2,ref_2] = ler(u(1),ref_1,tskip+delta,1,u(2));       % find the new reference state
        actuals(:,iEq) = [sigma_2 u(2)];
        hold on;
        scatter(u(2), std(getHeadways(ref_2(1:numCars))), 400, 'b.'); drawnow; %% plot the bifurcation diagram
 
    end
    save('origEqFree.mat', 'bif', 'actuals');

    
    %% functions
    
    %% function to zero for fsolve
    % u         - the current value of (sigma, v0) that we're trying to find
    %               with Newton's method
    % ref       - the most recent reference state
    % W         - the slope of the secant line for arc length continuation
    % newGuess  - the first guess on the secant line for arc length
    %               continuation
    %
    % RETURNS:
    % fw    - the functions F and w evaluated at these parameters
    function fw = FW(u,ref,W,newGuess)
        fw = zeros(2,1);
        fw(1) = F(ref,u(1),u(2));
        fw(2) = W(1)*(u(1)-newGuess(1)) + W(2)*(u(2) - newGuess(2));
    end

    %% lift, evolve, restrict
    % sigma - the current value of the std, used to seed the lifting
    % ref   - a reference state, used to seed the lifting
    % t     - the duration to evaluate the lifted parameters
    % p     -  parameter which will degrade the accuracy of the lifting
    %               operator.  Should be kept at p = 1 by default.
    % v0    - the optimal velocity parameter for this state
    % RETURNS: 
    % sigma     - the std. of the headways after restricting has occured
    % new_state - the final state of the evolution, which can be used as a
    %               future reference state
    function [sigma,new_state, sigma2, new_state2] = ler(sigma,ref,t,p,v0, tReference)
        lifted = lift(sigma, p, getHeadways(ref(1:numCars)),v0); %lift
         [~,evo] = ode45(@microsystem,[0 t],lifted, options,v0); % evolve
         if (nargin > 5) % evolve for tskip
             [~,evo2] = ode45(@microsystem,[0 tReference],evo(end,1:2*numCars)',options,v0);
            sigma2 = std(getHeadways(evo2(end,1:numCars)'));
         end
         evoCars = evo(end, 1:numCars)';
         sigma = std(getHeadways(evoCars)); % restrict
         if(nargout >= 2)
             new_state = evo(end,1:2*numCars)';
         end
         if(nargout == 4)
             new_state2 = evo2(end,1:2*numCars)';
         end
    end

    %% Finite Difference Quotient
    %  ref - reference state to base the lifting 
    %  sigma - the current value of the std. at which point to approximate
    %       the time derivative
    %  v0 - the velocity parameter for this state
    %  RETURNS:
    %  dif - the difference which approximates the time derivative
    function dif = F(ref, sigma,v0)
        [r0, ~, r1] = ler(sigma, ref, tskip, 1, v0, delta);
        dif = (r1-r0)/delta;
    end

    %% Jacobian for newton's method
    % ref - The previous reference state used to compute F.
    % sigma - the current value of sigma 
    % v0 - the velocity parameter for this state
    % w - the secant direction
    % J- The Jacobian, which will be given by
    % | F_sigma    F_v0  |
    % | w_sigma    w_vo  |
    function J = jacobian(ref, sigma, v0,w)
        J = zeros(2);
        unchanged = F(ref, sigma, v0);
        J(1,1) = (F(ref, sigma + delSigma, v0) - unchanged)/delSigma;
        J(1,2) = (F(ref, sigma,v0 + delv0) - unchanged)/delv0;
        J(2,:) = w';
    end
    
    %% lifting operator
    % s 	-   the target std. for this state
    % p     -   parameter which will degrade the accuracy of the lifting
    %           operator.  Should be kept at p = 1 by default.
    % hways -   headways of the cars at a reference state
    % v0    -   the goal velocity at this state
    % RETURNS: 
    % new   -   the lifted vector of cars and their positions. The cars
    %           will be reset to new positions with the first car centered 
    %            at 0, and with headways that will give the desired std. of
    %            s.  The cars' new velocities will be given directly by the
    %            optimal velocity function.
    function new = lift(s, p, hways, v0)
        del = p * s/std(hways) * (hways - mean(hways)) + mean(hways);
        new = zeros(length(hways)*2,1);
        for iLift=2:length(hways)
            new(iLift) = sum(del(1:iLift-1));
        end
        hways = getHeadways(new(1:numCars));
        new(numCars+1:end) = optimalVelocity(hways, v0);
    end
    
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