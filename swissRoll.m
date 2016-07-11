function swissRoll()

%% swill roll from Dsilva paper
% number of data points
N = 500;

% construct archemedian spiral
a = 1;
theta_vec = linspace(0, 4*pi, 100);
s = 0.5*a*(theta_vec.*sqrt(1+theta_vec.^2)+log(theta_vec+sqrt(1+theta_vec.^2)));
% height
h = 20;
%% generate data
% intialize random number generator
rng(321);
% find angles which correspond to uniform sampling along spiral
theta = interp1(s, theta_vec, rand(N, 1)*max(s));
% data uniformly distributed on swiss roll
z = h*rand(N,1);
x = a * cos(theta) .* theta;
y = a * sin(theta) .* theta;
% store all data
data = [x y z];

%{
% plot the data
figure;
scatter3(data(:,1), data(:,2), data(:,3));
title('Swiss roll data points');
%}

% calculate the distance matrix
D = zeros(N);
for r = 1:N
    for c = 1:N
        D(r,c) = norm(data(r,:)-data(c,:));
    end
end

% find the value of epsilon
eps = sqrt(median(pdist(data))/3);

% find the diffusion map
k = 5;
[vec, val] = diffusionMap(eps, D,k);
% 
 rk3 = linearFit(vec,3)
% rk4 = linearFit(vec,4)
% rk5 = linearFit(vec,5)
% 
% things = [rk3 rk4 rk5]

% plot the data colored by the eigen direction
figure;
scatter3(data(:,1), data(:,2), data(:,3), 200, vec(:,1),'.');
title('Data colored by first eigen-direction with h = 20');

figure;
scatter3(data(:,1), data(:,2), data(:,3), 200, vec(:,3),'.');
title('Data colored by fifth eigen-direction with h = 20');


    function plotEps()
        epsilon = 0:stepSize:maxEps;
        L = zeros(size(epsilon));
        for iEps = 1:length(epsilon)
            curEps = epsilon(iEps);
            L(iEps) = sum(sum(D<curEps,1),2);
        end
        
        figure;
        loglog(epsilon,L);
        xlabel('\epsilon','FontSize',24);
        ylabel('L(\epsilon)','FontSize',20);
    end

    function hways = getHeadways(v)
        futureCars = circshift(v,[-1,0]);
        hways = mod(futureCars - v, 60);
    end
end