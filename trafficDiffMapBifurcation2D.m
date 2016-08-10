function trafficDiffMapBifurcation2D()
h = 1.2;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 60;           % number of cars
stepSize = .0002;        % step size for the secant line approximation
tolerance = 10^(-14);    % tolerance for Newton's method
options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
foptions = optimset('TolFun',tolerance);              % fsolve options

%% load diffusion map data
load('bigDataMap.mat', 'evecs', 'evals', 'eps', 'allData');

% calculate which eigenvectors are most significant
%{
r = zeros(numEigvecs, 1);
r(1) = 1;
for j = 2:numEigvecs
    r(j) = linearFit(evecs,j);
end
%}
%{
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
%}

%% initialize secant continuation
steps = 100;                                % number of steps to take around the curve
bif = zeros(3,steps);                       % array to hold the bifurcation values

% initialize the first reference state
load('start0.884800.mat', 'trafficData', 'vel');
[trafficOutput, tan,t2] = findPeriodic(trafficData, evals, evecs, allData, eps, vel);
ref_2 = getHeadways(trafficOutput(1:60), len);
v0_base2 = vel;
leslieann_2 = diffMapRestrict(ref_2,evals,evecs,allData,eps);
start = leslieann_2;

% initialize the second reference state
load('start0.885000.mat', 'trafficData', 'vel');
[looped, ~, ~] = findPeriodic(trafficData, evals, evecs, allData, eps, vel, tan, start);
[~,~,t1] = findPeriodic(looped,evals,evecs,allData,eps,vel);
ref_1 = getHeadways(looped(1:60), len);
v0_base1 = vel;
leslieann_1 = diffMapRestrict(ref_1,evals,evecs,allData,eps);      %initial sigma values for secant line approximation
x0 = leslieann_1;
embed_1 = [0, t1, v0_base1];
embed_2 = [norm(leslieann_1-leslieann_2), t2, v0_base2];

psi = (leslieann_2 - leslieann_1)/norm(leslieann_2 - leslieann_1);

figure; hold on;
scatter(evecs(:,1), evecs(:,2), 'c.');
scatter(leslieann_1(1), leslieann_1(2), 'b*');
scatter(leslieann_2(1), leslieann_2(2), 'ro');
drawnow;

%% pseudo arc length continuation
for iEq=1:steps
    fprintf('Starting iteration %d of %d \n', iEq, steps);
    w = embed_2-embed_1;          % slope of the secant line
    newGuess = embed_2 + stepSize *(w/norm(w)); % first guess on the secant line
                                	
    %% alternate Newton's method using fsolve
    u = fsolve(@(u)FW(u,allData,w,newGuess,evecs,evals,eps), newGuess,foptions)
    
    bif(:,iEq) = u;                                            % save the new solution
    
    %% reset the values for the arc length continuation
    leslieann_2 = u(1)*psi + x0;                     % find the new reference state
    scatter(leslieann_2(1),leslieann_2(2),'m.'); drawnow;
    
    embed_1 = embed_2;
    embed_2 = u;
end
hold off;

%% plot the bifurcation diagram
figure;
scatter3(bif(1,:),bif(2,:),bif(3,:),'*');
xlabel('\alpha');
ylabel('T');
zlabel('v_0');

figure;
scatter(bif(3,:), bif(1,:), 300, 'b.');
xlabel('v_0');
ylabel('\alpha');

%% Jacobian for newton's method
% ref - The previous reference state used to compute F.
% sigma - the current value of sigma
% v0 - the velocity parameter for this state
% w - the secant direction
% J- The Jacobian, which will be given by
% | F_sigma    F_v0  |
% | w_sigma    w_vo  |
    function J = jacobian(ref, embed, v0,w,eigvecs,eigvals,lereps)
        J = zeros(3);
        delAlpha = 10^(-7);
        delT = 10^(-3);
        delv0 = 10^(-7);
        x1 = embed(1)*psi + x0;
        xDelta = (embed(1) + delAlpha)*psi + x0;
        unchanged = F(ref, x1, v0,eigvecs,eigvals,lereps, embed(2));
        J(1:2,1) = (F(ref, xDelta, v0,eigvecs,eigvals,lereps, embed(2)) - unchanged)/delAlpha;
        J(1:2,2) = (F(ref, x1,v0,eigvecs,eigvals,lereps,embed(2)+delT) - unchanged)/delT;
        J(1:2,3) = (F(ref, x1,v0 + delv0,eigvecs,eigvals,lereps, embed(2)) - unchanged)/delv0;
        J(3,:) = w';
    end

%% function to zero for fsolve
% u         - the current value of (alpha,T,v0)
% ref       - the most recent reference state
% W         - the slope of the secant line for arc length continuation
% newGuess  - the first guess on the secant line for arc length
%               continuation
%
% RETURNS:
% fw    - the functions F and w evaluated at these parameters
    function fw = FW(u,ref,W,newGuess,evecs,evals,lereps)
        x = u(1)*psi + x0;
        scale = [1 ; 1 ; 1];
        fw = zeros(3,1);
        fw(1:2) = F(ref,x(1:2),u(end),evecs,evals,lereps,u(2)) .* scale(1:2);
        fw(end) = dot(W,u - newGuess) * scale(end);
        fprintf('|F| = %d\n\n',norm(fw));
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
    function [sigma, sigma2] = ler(newval,orig,t,v0,eigvecs,eigvals,lereps,tReference)
        fprintf('Lifting to:  %f \n', newval(1));
        fprintf('\t and %f \n', newval(2));
        lifted = smartLift2d(newval, eigvecs, eigvals, lereps,v0, orig);
        [~,evo] = ode45(@microsystem,[0 t],lifted, options,[v0 len h]);        
        evoCars = getHeadways(evo(end,1:numCars)',len);
        sigma = diffMapRestrict(evoCars, eigvals, eigvecs, orig, lereps);
        if (nargin > 7)
            [~,evo2] = ode45(@microsystem,[0 tReference],evo(end,:)', options,[v0 len h]);
            evo2Cars = getHeadways(evo2(end,1:numCars)',len);       
            sigma2 = diffMapRestrict(evo2Cars,eigvals,eigvecs, orig, lereps);
        end
    end

%% Finite Difference Quotient
%  ref - reference state to base the lifting
%  sigma - the current value of the std. at which point to approximate
%       the time derivative
%  v0 - the velocity parameter for this state
%  RETURNS:
%  dif - the difference which approximates the time derivative
    function dif = F(ref, sigma,v0,eigvecs,eigvals,lereps,t)
        delta = 28*t;
        tskip = 4*t;
        [r0, r1] = ler(sigma, ref, tskip, v0, eigvecs,eigvals,lereps, delta);
        dif = (r1-r0)/delta;
    end

end