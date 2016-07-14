%% The function runDiffMap takes in data from traffic simulations,
% builds a diffusion map, compares the eigen directions given by the
% diffusion map eigenvectors, and plots relevent eigenvectors scaled by
% corresponding eigenvalues
function runDiffMap()
%jacobOutput and jacobOutput2 only contain data for fully developed
%simulations, with a constant v0 of 0.885 and final time of t = 70000
% and sin amplitude between 0.55 and 1.
load('jacobOutput','joutput885');
load('jacobOutput2','j2output885');
load('jacobOutput3', 'j3output885');
load('jacobOutput7', 'j7output885');

%jacobOutput 3 contains data for constant v0 of 0.885, but with a final
%time uniformly distributed between t = 35000 and 55000

%jacobOutput7 contains data for constant v0 of 0.885, but with a final
%time uniformly distributed between t = 35000 and 45000 and magnitude
%between 0.75 and 1

% build the data from the simulation results
allTime = [joutput885 j2output885 j3output885 j7output885];
allTime = getHeadways(allTime(1:60,:)); % find the headways from positions
[~, max1] = max(allTime(1:60,:),[],1);  % locate the max headway for each data point

% align all of the headways with the max in the front
for c = 1:length(max1)
    allTime(1:60,c) = circshift(allTime(1:60,c), [-max1(c)+1,0]);
end

newData = allTime(:,end-49:end);
newstdev = std(newData, 0, 1);
allTime = allTime(:,1:end-50);
stdev = std(allTime, 0, 1);             % find the standard deviations

% calcuate the pairwise distances between data points
D = zeros(length(allTime));
for r = 1:length(allTime)
    for c = 1:length(allTime)
        D(r,c) = norm(allTime(:,r)-allTime(:,c));
    end
end

eps = median(D(:)); % choose epsilon for the kernel based on the pairwise distances

k=1;      % number of eigenvectors to calculate from the diffusion map
t = 3;      % number of timesteps to evolve the diffusion map
[vec,val] = diffusionMap(eps,D,k);          % calculate the diffusion map

dvec = zeros(k,50);
for iPt = 1:50
    dvec(:,iPt) = diffMapRestrict(newData(:,iPt),val,vec,allTime,eps);
end

diffMapRestrict(diffMapLift(vec(501,1),vec(:,1),allTime),val,vec,allTime,eps)
vec(501,1)
% calculate how unique each eigen direction is
%{
r = zeros(k, 1);
r(1) = 1;
for j = 2:k
    r(j) = linearFit(vec,j);
end
% display the eigen direction computation results
fprintf('The r_k values are: \n');
disp(r);
%}

% evolve the diffusion map t times
val = val^t;

%Change this to see the relationship between the eigenfunction at eigDisp
%and the standard deviation.

eigDisp = 1;
stry = sprintf('Eigenvector %i',eigDisp);
figure;
hold on;
scatter(stdev, vec(:,eigDisp),'b.')
scatter(newstdev, dvec(eigDisp, :),500,'r.')
xlabel('\sigma');
ylabel(stry);


% Plot relationship between first three eigenvectors
figure;
hold on;
scatter(vec(:,1), vec(:,2),100,'b.');
scatter(dvec(1,:),dvec(2,:),100,'r.');
hold off;
xlabel('Eigenvector 1');
ylabel('Eigenvector 2');


    % plotEps creates a log-log plot of the number of data points that are
    % less than epsilon vs. epsilon
    function plotEps(distances)
        % create the values of epsilon to test from 0 to maxEps by stepSize
        stepSize = .0001;
        maxEps = 1;
        epsilon = 0:stepSize:maxEps;
        % count the number of points that are less than each epsilon
        L = zeros(size(epsilon));
        for iEps = 1:length(epsilon)
            curEps = epsilon(iEps);
            L(iEps) = sum(sum(distances<curEps,1),2);
        end
        % plot the log-log plot of L(epsilon) vs epsilon
        figure;
        loglog(epsilon,L);
        xlabel('\epsilon','FontSize',24);
        ylabel('L(\epsilon)','FontSize',20);
    end

    % getHeadways returns headways between 60 cars given their positions
    function hways = getHeadways(v)
        futureCars = circshift(v,[-1,0]);
        hways = mod(futureCars - v, 60);
    end
end