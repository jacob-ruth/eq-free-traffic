function trafficDiffMapBifurcation2D()
h = 1.2;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 60;           % number of cars
tskip = 300;            % times for evolving
delta = 500;
stepSize = .0002;        % step size for the secant line approximation
tolerance = 10^(-12);    % tolerance for Newton's method
options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
foptions = optimset('TolFun',tolerance);              % fsolve options

%% load diffusion map data
load('bigDataMap.mat', 'evecs', 'evals', 'eps', 'allData');

[~, max1] = max(allData,[],1);  % locate the max headway for each data point
% plot eigenvector 1 vs eigenvector 2, colored by max headway location
figure;
scatter(evecs(:,1), evecs(:,2), 100, max1,'.');
h = colorbar;
xlabel(h, 'Wave Position', 'fontsize', 15)
colormap(jet);
xlabel('\Phi_1', 'FontSize',15);
ylabel('\Phi_2', 'FontSize',15);
title('\Phi_1 vs. \Phi_2 Colored by Locations of the Max Headways','FontSize',15);

% plot eigenvector 1 vs eigenvector 2, colored by standard deviation
figure;
scatter(evecs(:,1), evecs(:,2), 100,  std(allData),'.');
colorbar;
h = colorbar;
xlabel(h, '\sigma', 'fontsize', 15)
colormap(jet);
xlabel('\Phi_1', 'FontSize',15);
ylabel('\Phi_2', 'FontSize',15);
title('\Phi_1 vs. \Phi_2 Colored by Standard Deviation of the Headways','FontSize',15);

%% initialize secant continuation
steps = 100;                                % number of steps to take around the curve
bif = zeros(3,steps);                       % array to hold the bifurcation values

% initialize the first reference state
load('start0.884800.mat', 'trafficData', 'vel');
[trafficOutput, tang] = findPeriodic(trafficData, evals, evecs, allData, eps, vel);
ref_2 = getHeadways(trafficOutput(1:60), len);
v0_base2 = vel;
embed_2 = diffMapRestrict(ref_2,evals,evecs,allData,eps);
start = embed_2;

% initialize the second reference state
load('start0.885000.mat', 'trafficData', 'vel');
looped = findPeriodic(trafficData, evals, evecs, allData, eps, vel, tang, start);
ref_1 = getHeadways(looped(1:60), len);
v0_base1 = vel;
embed_1 = diffMapRestrict(ref_1,evals,evecs,allData,eps);      %initial sigma values for secant line approximation

figure; hold on;
scatter(evecs(:,1), evecs(:,2), 'c.');
scatter(embed_1(1), embed_1(2), 'b*');
scatter(embed_2(1), embed_2(2), 'ro');

%% pseudo arc length continuation
for iEq=1:steps
    fprintf('Starting iteration %d of %d \n', iEq, steps);
    w = [embed_2 - embed_1 ; v0_base2 - v0_base1];          % slope of the secant line
    newGuess = [embed_2; v0_base2] + stepSize *(w/norm(w)); % first guess on the secant line
                                	
    %% alternate Newton's method using fsolve
    u = fsolve(@(u)FW(u,allData,w,newGuess,evecs,evals,eps), newGuess,foptions);
    
    bif(:,iEq) = u;                                            % save the new solution
    scatter(u(1), u(2), 'k.'); drawnow;
    
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
ylabel('\Phi_2');
zlabel('v_0');

rad = sqrt(bif(1,:).^2 + bif(2,:).^2);
figure;
scatter(bif(3,:), rad, 'b.');
xlabel('v_0');
ylabel('\sigma');

%% function to zero for fsolve
% u         - the current value of (phi_1, v0) that we're trying to find
%               with Newton's method
% ref       - the most recent reference state
% W         - the slope of the secant line for arc length continuation
% newGuess  - the first guess on the secant line for arc length
%               continuation
%
% RETURNS:
% fw    - the functions F and w evaluated at these parameters
    function fw = FW(u,allData,W,newGuess,evecs,evals,lereps)
        fw = zeros(3,1);
        fw(1:2) = F(allData,u(1:2),u(end),evecs,evals,lereps);
        fw(end) = W'*(u - newGuess);
        fprintf('F = %d \n\n',norm(fw));
    end

%% lift, evolve, restrict
% newval     - the current value of the macrovariables,
%                 used to seed the lifting
% orig       - the data set, for use in lifting
% t          - the duration to evaluate the lifted parameters
% v0         - the optimal velocity parameter for this state
% eigvecs    - eigenvectors from diffusion map
% eigvals    - eigenvalues from diffusion map
% lereps     - epsilon used in diffusionmap
% tReference - a second state to lift to, if desired
% RETURNS:
% sigma      - the macro-state of the headways after evolving for t
% sigma2     - the macro-state of the headways, after evolving for t +
%               tReference
    function [sigma, sigma2] = ler(newval,orig,t,v0,eigvecs,eigvals,lereps,tReference)
        fprintf('Lifting to: %f \n', newval(1));
        fprintf('\t and %f \n', newval(2));
        lifted = smartLift2d(newval, eigvecs, eigvals, lereps,v0, orig);
        [~,evo] = ode45(@microsystem,[0 t],lifted, options,[v0 len h]);
        evo = findPeriodic(evo(end,:)', eigvals, eigvecs, orig, lereps, v0, tang, start);        
        evoCars = getHeadways(evo(1:numCars),len);
        sigma = diffMapRestrict(evoCars, eigvals, eigvecs, orig, lereps);
        if (nargin > 7)
            [~,evo2] = ode45(@microsystem,[0 tReference],evo, options,[v0 len h]);
            evo2 = findPeriodic(evo2(end,:)', eigvals, eigvecs, orig, lereps, v0, tang, start);
            evo2Cars = getHeadways(evo2(1:numCars),len);       
            sigma2 = diffMapRestrict(evo2Cars,eigvals,eigvecs, orig, lereps);
        end
    end

%% Finite Difference Quotient
% allData    - data set for use in lifting
% sigma      - the current value of the macrovariable at which point to approximate
%               the time derivative
% v0         - the velocity parameter for this state
% eigvecs    - eigenvectors from diffusion map
% eigvals    - eigenvalues from diffusion map
% lereps     - epsilon used in diffusionmap

%  RETURNS:
%  dif - the difference which approximates the time derivative
    function dif = F(allData, sigma,v0,eigvecs,eigvals,lereps)
        [r0, r1] = ler(sigma, allData, tskip, v0, eigvecs,eigvals,lereps, delta);
        dif = (r1-r0)/delta;
    end

end