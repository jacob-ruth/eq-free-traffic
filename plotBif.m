%% plot the bifurcation diagram
% load the diffusion map bifurcation data
load('diffusionMapBif.mat', 'x', 'y');

% load the sigma bifurcation data
load('eqFreeBif.mat', 'bif');

% plot the two bifurcation diagrams
figure; hold on;
plot(x, y, 'r.', 'markersize', 12);
plot(bif(2,30:end), bif(1,30:end), 'b.', 'markersize', 12);
xlabel('v_0', 'fontsize', 16);
ylabel('\sigma', 'fontsize', 20);
title('Bifurcation Diagram', 'fontsize', 20);
ylim([0 0.2]);
legend('Diffusion Map Embedding', 'Original Equation-Free Approach');
hold off;