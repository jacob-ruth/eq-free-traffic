 function trafficDiffMapNewtonContinuation2D()
clear workspace;

h = 2.4;                % optimal velocity parameter
len = 60;               % length of the ring road
numCars = 30;           % number of cars
lsqOptions = optimset('Display','iter'); %lsqnonlin options
options2 = optimoptions('lsqnonlin');
options2.OptimalityTolerance = 1e-10;
options2.FunctionTolerance = 1e-10;
options2.Display = 'iter';

%% load diffusion map data
%load('bigDataMap.mat', 'evecs', 'evals', 'eps', 'allData');

%load('newDiffMap.mat', 'hWays', 'data', 'evecs', 'evals', 'eps')
load('30DiffMap.mat', 'hways', 'eps', 'evecs', 'evals');
allData=hways;

%{
[~, max1] = max(allData,[],1);  % locate the max headway for each data point
% plot eigenvector 1 vs eigenvector 2, colored by max headway location
figure(1);
scatter(evecs(:,1), evecs(:,2), 100, max1,'.'); hold on;
color = colorbar;
xlabel(color, 'Wave Position', 'fontsize', 16)
colormap(jet);
xlabel('\psi_1', 'FontSize',18);
ylabel('\psi_2', 'FontSize',18);
title('\psi_1 vs. \psi_2 Colored by Locations of the Max Headways','FontSize',16);
drawnow;


% plot eigenvector 1 vs eigenvector 2, colored by standard deviation
figure;
scatter(evecs(:,1), evecs(:,2), 100,  std(allData),'.');
colorbar;
color = colorbar;
xlabel(color, '\sigma', 'fontsize', 18)
colormap(jet);
xlabel('\psi_1', 'FontSize',18);
ylabel('\psi_2', 'FontSize',18);
title('\psi_1 vs. \psi_2 Colored by Standard Deviation of the Headways','FontSize',16);
%}


addKernelData(allData, eps); % save kernel data to avoid computing every time

%{
% initialize the first reference state
load('start0.884800.mat', 'trafficData', 'vel');
[trafficOutput, tang] = findPeriodic_bjorn(trafficData, evals, evecs, allData, eps, vel); % finds the evolved state that has gone around the ring road exactly once
ref_2 = getHeadways(trafficOutput(1:60), len);
v0_base2 = vel;
embed_2 = diffMapRestrictAlt(ref_2,evals,evecs,allData,eps,1);
% each triple is (sigma, T, v0)
p2 = norm(embed_2);
T2 = travelingWavePeriod_bjorn(ref_2);
coord2 = [p2, T2, v0_base2];
rayAngle = atan2(embed_2(2),embed_2(1));
start = embed_2;

% initialize the second reference state
load('start0.885000.mat', 'trafficData', 'vel');
looped = findPeriodic_bjorn(trafficData, evals, evecs, allData, eps, vel, tang, start);
ref_1 = getHeadways(looped(1:60), len);
v0_base1 = vel;
embed_1 = diffMapRestrictAlt(ref_1,evals,evecs,allData,eps,1);
p1 = norm(embed_1);
T1 = travelingWavePeriod_bjorn(ref_1);
coord1 = [p1, T1, v0_base1];
%}


load('microBif.mat', 'bif', 'vel',  'period', 'n');

start = 110;
change = 1;
v0_base2 = vel(start + change);
v0_base1 = vel(start);
ref_2 = bif(1:numCars,start + change);
ref_1 = bif(1:numCars,start);
%ref_2 = alignMax(bif(1:numCars,start + change),10,0);
%ref_1 = alignMax(bif(1:numCars,start),10,0);
T2 = -numCars/bif(numCars + 1, start + change);
T1 = -numCars/bif(numCars + 1, start);
embed_2 = diffMapRestrictAlt(ref_2,evals,evecs,allData,eps,1);
embed_1 = diffMapRestrictAlt(ref_1,evals,evecs,allData,eps,1);
p2 = norm(embed_2);
p1 = norm(embed_1);
%}
%{
input = embed_1;
output = ler(input, allData, T1, v0_base1, evecs, evals, eps);
output2 = ler(input, allData, period(start), v0_base1, evecs, evals, eps);
%[T,looped] = findPeriodic_bjorn(ref_1, [hwayToPos(ref_1) ; optimalVelocity(h,ref_1, v0_base1)], evals, evecs, allData, eps, v0_base1);
%output3= ler(input, allData,T(end)-T(end-1), v0_base1, evecs, evals, eps);
%hw = getHeadways(looped(1:numCars),len);
%embed = diffMapRestrictAlt(hw,evals,evecs,allData,eps,1);
%output4 = ler(embed, allData,T(end)-T(end-1), v0_base1, evecs, evals, eps);
figure; hold on;
scatter(evecs(:,1), evecs(:,2), 50, 'k.');
scatter(input(1), input(2), 400, 'r.');
scatter(output(1), output(2), 400, 'b.');
scatter(output2(1), output2(2), 400, 'g.');
%scatter(output3(1), output3(2), 400, 'y.');
%scatter(embed(1), embed(2), 200, 'r*');
%scatter(output4(1), output4(2), 200, 'b*');
%}
%{
load('profile2.mat', 'allTime', 'v0');
ref_2 = getHeadways(allTime(1:numCars)', len);
cars2 = allTime;
v0_base2 = v0;
load('profile1.mat', 'allTime', 'v0');
ref_1 = getHeadways(allTime(1:numCars)', len);
v0_base1 = v0;
cars1 = allTime;
T2 = findPeriodic_bjorn(ref_2, cars2, evals, evecs, allData, eps, v0_base2);
T1 = findPeriodic_bjorn(ref_1, cars1, evals, evecs, allData, eps, v0_base1);
%}

coord2 = [p2, T2, v0_base2];
coord1 = [p1, T1, v0_base1];
rayAngle = atan2(embed_2(2), embed_2(1));

figure;
subplot(1,2,1);hold on;
scatter(vel(100:end), n(100:end), 50, 'r.');
xlabel('v_0');
ylabel('\rho');
scatter(v0_base1, p1, 400,'k.'); drawnow;
scatter(v0_base2, p2, 400,'k.'); drawnow;

subplot(1,2,2); hold on;
scatter(n(100:end), -numCars./bif(numCars + 1, 100:end), 50, 'r.');
scatter(p1,T1, 400, 'k.'); drawnow;
scatter(p2, T2, 400, 'k.'); drawnow;
xlabel('\rho');
ylabel('T');

%% initialize secant continuation
steps = 500;                                % number of steps to take around the curve
bif = zeros(3,steps);                      % array to hold the bifurcation values
sigma = zeros(1,steps);
guesses = zeros(3,steps);
stepSize = .003;        % step size for the secant line approximation

%% pseudo arc length continuation
for iEq=1:steps
    fprintf('Starting iteration %d of %d \n', iEq, steps);
    w = coord2 - coord1;
    scaledW = w ./ coord2;
    newGuess = coord2 + stepSize *(w/norm(w)); % first guess on the secant line
    guesses(:,iEq) = newGuess;
   
    scaledNewGuess = newGuess ./ coord2;

    %% alternate Newton's method using lsqnonlin
    [u, ~, ~, exitFlag] = lsqnonlin(@(u)periodicDistance(u,allData,evecs,evals,eps,coord2,rayAngle,scaledW,...
        scaledNewGuess), newGuess,[],[],lsqOptions);
    
    
    bif(:,iEq) = u';                                            % save the new solution
    
    %% reset the values for the arc length continuation
    coord1 = coord2;
    coord2 = u;
    
    sigma(iEq) = sig;
    
    %fprintf('New solution:    %d, %d, %d %d \n', u(1), u(2), u(3), sigma(iEq));
    fprintf('Exit flag: %d \n', exitFlag);
    save('newtonContinuation7.mat', 'bif', 'guesses', 'rayAngle', 'sigma');
    
    subplot(1,2,1); hold on;
    scatter(u(3), u(1), 400, 'b.'); drawnow;
    subplot(1,2,2); hold on;
    scatter(u(1), u(2), 400, 'b.'); drawnow;
end

disp('done');

%% function to zero
% PARAMETERS:
% u         - the pair of radius and period (sigma,T,v0) that we're trying to find with
% allData
% evecs
% evals
% lereps
% coord2
% ray
% scaledW
% scaledNewGuess
% RETURNS:
% distance  - the distance between the input point and the point determined by 
%               lifting, evolving for T time, and restricting
% onplane   - whether this is on the plane orthogonal to w at the starting guess
    function out = periodicDistance(u,allData,evecs,evals,lereps, coord2, ray, scaledW, scaledNewGuess)
        inputPoint = [u(1)*cos(ray) ; u(1)*sin(ray)];
        [finalPoint, sig] = ler(inputPoint,allData,25*u(2),u(3),evecs,evals,lereps);
        
        angle = atan2(finalPoint(2), finalPoint(1));
        r = norm(finalPoint);
        distance = [((u(1)-r)/10^(-3))^2 ; (angle-ray)^2 ];

        scaledU = u ./ coord2;
        onPlane = dot(scaledW, scaledU - scaledNewGuess);
        out = [distance' onPlane];
        
        %{
        figure; hold on;
        scatter(evecs(:,1), evecs(:,2), 50, 'k.');
        scatter(inputPoint(1), inputPoint(2), 400, 'r.');
        scatter(finalPoint(1), finalPoint(2), 400, 'b.');
        pause;
        %}
    end


%% lift, evolve, restrict
% PARAMETERS:
% newval     - the current value of the macrovariables,
%                 used to seed the lifting
% orig       - the data set, for use in lifting
% t          - the duration to evaluate the lifted parameters
% v0         - the optimal velocity parameter for this state
% eigvecs    - eigenvectors from diffusion map
% eigvals    - eigenvalues from diffusion map
% lereps     - epsilon used in diffusionmap
% RETURNS:
% sigma      - the macro-state of the headways after evolving for t
    function [coord, sigma, evo] = ler(newval,orig,t,v0,eigvecs,eigvals,lereps)
        options = odeset('AbsTol',10^-8,'RelTol',10^-8); % ODE 45 options
        lifted = newLift(newval, eigvecs, eigvals, lereps,v0, orig);
        [~,evo] = ode45(@microsystem,[0 t],lifted, options,[v0 len h]);
        evo = evo(end,:)';
        evoCars = getHeadways(evo(1:numCars),len);
        sigma = std(evoCars);
        coord = diffMapRestrictAlt(evoCars,eigvals,eigvecs, orig, lereps, 1);
    end

end