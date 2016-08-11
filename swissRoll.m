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
    h = 40;                                         % height
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
    eps = median(D(:))/8;

    % find the diffusion map
    k = 5;                               % number of eigenvectors to return
    [vec, val] = diffusionMap(eps, D, k);
    
    specialPoint = data(1000,:);
    eigCoor = vec(1000,:)
    newData = [data(1:999,:) ; data(1001:end,:)];
    [newvec, newval] = diffusionMap(eps, squareform(pdist(newData)),k);
    newEigCoor = diffMapRestrict(specialPoint',newval,newvec,newData',eps)

%{
% compute how unique the eigen directions given by the vectors are
 r = zeros(k,1);
r(1) = 1;                    	% assume the first direction is important
for i=2:k
    r(i) = linearFit(vec, i);
end
% display the values of r
fprintf('The values of r_k are \n');
disp(r);
%}

    % plot the data colored by the eigen directions
    figure;
    scatter3(data(:,1), data(:,2), data(:,3), 500, vec(:,1),'.');
    title('Data Colored by First Eigen-Direction','FontSize',20);
    colormap(jet);
    h = colorbar;
    xlabel(h, '\Phi_1', 'fontsize', 15)

    figure;
    scatter3(data(:,1), data(:,2), data(:,3), 500, vec(:,2),'.');
    title('Data Colored by Second Eigen-Direction','FontSize',20);
    colormap(jet);
    h = colorbar;
    xlabel(h, '\Phi_2', 'fontsize', 15)

end
