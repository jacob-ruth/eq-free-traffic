%% Function swissRoll samples data points from a swiss roll, constructs a 
% diffusion map, compares the eigenvector directions, and plots the data 
% colored by the eigenvectors
function swissRoll()
    %% swiss roll from Dsilva paper
    N = 1500;                                       % number of data points
    % construct archemedian spiral 
    a = 1;
    theta_vec = linspace(0, 4*pi, 100);
    s = 0.5*a*(theta_vec.*sqrt(1+theta_vec.^2)+log(theta_vec+sqrt(1+theta_vec.^2)));
    h = 20;                                         % height
    %% generate data
    rng(321);                             % intialize random number generator
    % find angles which correspond to uniform sampling along spiral
    theta = interp1(s, theta_vec, rand(N, 1)*max(s));
    % data uniformly distributed on swiss roll
    z = h*rand(N,1);
    x = a * cos(theta) .* theta;
    y = a * sin(theta) .* theta;
    data = [x y z];                             % store all data

    % calculate the distance matrix
    D = squareform(pdist(data));

    % find the value of epsilon - sqrt(5) was used by Dsliva code
    eps = sqrt(5);

    % find the diffusion map
    k = 5;                               % number of eigenvectors to return
    [vec, ~] = diffusionMap(eps, D, k);

    % compute how unique the eigen directions given by the vectors are
    r = zeros(k,1);
    r(1) = 1;                    	% assume the first direction is important
    for i=2:k
        r(i) = linearFit(vec, i);
    end
    % display the values of r
    fprintf('The values of r_k are \n');
    disp(r);

    % plot the data colored by the eigen directions
    figure;
    scatter3(data(:,1), data(:,2), data(:,3), 200, vec(:,1),'.');
    title('Data colored by first eigen-direction');

    figure;
    scatter3(data(:,1), data(:,2), data(:,3), 200, vec(:,2),'.');
    title('Data colored by second eigen-direction');

    figure;
    scatter3(data(:,1), data(:,2), data(:,3), 200, vec(:,3),'.');
    title('Data colored by third eigen-direction');

    figure;
    scatter3(data(:,1), data(:,2), data(:,3), 200, vec(:,4),'.');
    title('Data colored by fourth eigen-direction');

    figure;
    scatter3(data(:,1), data(:,2), data(:,3), 200, vec(:,5),'.');
    title('Data colored by fifth eigen-direction');

end