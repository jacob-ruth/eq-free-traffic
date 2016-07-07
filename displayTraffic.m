function displayTraffic(u, L)
    [numVars,N] = size(u);
    carPositions = u(1:numVars/2, :);
    carPositions = rem(carPositions, L);
    stepSize = 5;
    close all; clc
    figure;

    % specify the plotting colors
    rbg = zeros(numVars/2,3);
    dx = 0.8;
    for i=1:numVars/2
        f = i*2/numVars;
        g = (6-2*dx)*f+dx;
        rgb(i,1) = max(0,(3-abs(g-4)-abs(g-5))/2);
        rgb(i,2) = max(0,(4-abs(g-2)-abs(g-4))/2); 
        rgb(i,3) = max(0,(3-abs(g-1)-abs(g-2))/2);
    end
    
    % Preallocate movie structure.
    mov(1:floor(N/stepSize)) = struct('cdata', [],'colormap', []);
    %set(gca,'nextplot','replacechildren');

    % plot the points
    for i = 1:stepSize:N
        angle = carPositions(:, i)*(2*pi)/L;
        x = cos(angle);
        y = sin(angle);
        clf;
        scatter(x,y,[],rgb);
        shg;
        pause(.01)
        mov(i) = getframe(gcf);
    end
    
    % create avi file
    movie2avi(mov, 'traffic.avi', 'compression', 'None');
end

