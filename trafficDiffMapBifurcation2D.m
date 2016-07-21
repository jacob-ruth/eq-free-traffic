function trafficDiffMapBifurcation2D()
h = 1.2;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 60;           % number of cars
tskip = 300;            % times for evolving
delta = 500;
stepSize = .0001;        % step size for the secant line approximation
delSigma = 0.00001;     % delta sigma used for finite difference of F
delv0 = 0.00001;        % delta v0 used for finite difference of F
tolerance = 10^(-5);    % tolerance for Newton's method
options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
foptions = optimset('TolFun',1e-8);              % fsolve options

%% load diffusion map data
load('jacobOutput','joutput885');
load('jacobOutput2','j2output885');
load('jacobOutput3', 'j3output885');
load('jacobOutput7', 'j7output885');
load('saveData','trafficOutput');
loadedData = [joutput885 j2output885 j3output885 j7output885 trafficOutput];
load('885low','trafficOutput');
loadedData = [loadedData trafficOutput];

allData = getHeadways(loadedData(1:numCars, :));        % headway data for the diffusion map

numEigvecs = 2;                                         % number of eigenvectors to return
[evecs, evals, eps] = runDiffMap(allData,numEigvecs);   % run the diffusion map

%{
% calculate which eigenvectors are most significant
r = zeros(numEigvecs, 1);
r(1) = 1;
for j = 2:numEigvecs
    r(j) = linearFit(evecs,j);
end
r
%}

% plot eigenvector 1 vs eigenvector 2, colored by max headway location
figure;
scatter(evecs(:,1), evecs(:,2), 10, max1);
colorbar;
xlabel('\Phi_1');
ylabel('\Phi_2');
title('Colored by max');

% plot eigenvector 1 vs eigenvector 2, colored by standard deviation
figure;
scatter(evecs(:,1), evecs(:,2), 10, std(allData));
colorbar;
xlabel('\Phi_1');
ylabel('\Phi_2');
title('Colored by standard deviation');

%% initialize secant continuation
steps = 50;                                % number of steps to take around the curve
bif = zeros(3,steps);                       % array to hold the bifurcation values

% initialize the first reference state
load('save884','trafficOutput','v0');
ref_2 = getHeadways(trafficOutput(1:60));
v0_base2 = v0;

% initialize the second reference state
load('885ref','trafficOutput','v0');
ref_1 = getHeadways(trafficOutput(1:60));
v0_base1 = v0;

embed_1 = diffMapRestrict(ref_1,evals,evecs,allData,eps);      %initial sigma values for secant line approximation
embed_2 = diffMapRestrict(ref_2,evals,evecs,allData,eps);

%% pseudo arc length continuation
for iEq=1:steps
    fprintf('Starting iteration %d of %d \n', iEq, steps);
    w = [embed_2 - embed_1 ; v0_base2 - v0_base1];          % slope of the secant line
    newGuess = [embed_2; v0_base2] + stepSize *(w/norm(w)); % first guess on the secant line
    
    %% initialize Newton's method
    u = newGuess;
    first = true;                       % mimic a do-while loop
    k=1;                            	% Newton's method counter
    
    %% Newton and that other guy's method
    while(first ||(norm(invD*[f;neww])>tolerance && k < 20))
        first = false;
        fprintf('\t Newton iteration: %d \n', k);
        f = F(allData, u(1),u(2), evecs, evals, eps);                                % calculcate the function to zero
        neww = w(1)*(u(1)-newGuess(1)) + w(2)*(u(2) - newGuess(2));
        Df = jacobian(allData, u(1), u(2), w, evecs, evals, eps);                    % find the jacobian
        invD = Df^(-1);
        u = u - invD*[f;neww]                                  % perform the Newton step
        k = k + 1;
        end
    
        %% alternate Newton's method using fsolve
    u = fsolve(@(u)FW(u,allData,w,newGuess,evecs,evals,eps), newGuess,foptions)
    
    bif(:,iEq) = u;                                            % save the new solution
    
    
    %% reset the values for the arc length continuation
    embed_1 = embed_2;
    v0_base1 = v0_base2;
    v0_base2 = u(end);
    embed_2 = u(1:2);                     % find the new reference state
end
hold off;

%% plot the bifurcation diagram
figure;
scatter3(bif(1,:),bif(2,:),bif(3,:),'*');
xlabel('\Phi_1');
ylabel('Phi_2');
zlabel('v0');

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
    function fw = FW(u,ref,W,newGuess,evecs,evals,lereps)
        fw = zeros(2,1);
        fw(1) = F(ref,u(1),u(2),evecs,evals,lereps);
        fw(2) = W'*(u - newGuess);
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
    function [sigma,new_state, sigma2, new_state2] = ler(newval,orig,t,v0,eigvecs,eigvals,lereps,tReference)
        liftedHways = diffMapLift(newval, eigvecs,eigvals,lereps, orig);
        liftedPosns = cumsum(liftedHways);
        liftedVel = optimalVelocity(liftedHways,v0);
        lifted = [liftedPosns ; liftedVel];
        [~,evo] = ode45(@microsystem,[0 t],lifted, options,v0);
        if (nargin > 7)
            [~,evo2] = ode45(@microsystem,[0 tReference],evo(end,1:2*numCars)',options,v0);
            evo2Cars = evo2(end, 1:numCars)';       
            evo2Cars = shiftMax(getHeadways(evo2Cars));
            sigma2 = diffMapRestrict(evo2Cars,eigvals,eigvecs, orig, lereps);
        end
        evoCars = evo(end, 1:numCars)';
        evoCars = shiftMax(getHeadways(evoCars));
        sigma = diffMapRestrict(evoCars, eigvals, eigvecs, orig, lereps);
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
    function dif = F(ref, sigma,v0,eigvecs,eigvals,lereps)
        [r0, ~, r1] = ler(sigma, ref, tskip, v0, eigvecs,eigvals,lereps, delta);
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
    function J = jacobian(ref, sigma, v0,w,eigvecs,eigvals,lereps)
        J = zeros(2);
        unchanged = F(ref, sigma, v0,eigvecs,eigvals,lereps);
        J(1,1) = (F(ref, sigma + delSigma, v0,eigvecs,eigvals,lereps) - unchanged)/delSigma;
        J(1,2) = (F(ref, sigma,v0 + delv0,eigvecs,eigvals,lereps) - unchanged)/delv0;
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
            new(iLift) = mod(sum(del(1:iLift-1)),len);
        end
        hways = getHeadways(new(1:length(hways)));
        for iLift=1:length(hways)
            new(iLift+length(hways)) = optimalVelocity(hways(iLift), v0);
        end
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