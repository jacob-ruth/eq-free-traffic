%% restriction operator based on diffusion map
% newData   - point to get restriction values for
% evals     - square matrix of eigenvalues
% evecs     - matrix with eigenvectors in columns
% origData  - matrix of data in rows used to compute eigenvalues/vectors
% eps       - length scale for kernel
% 
% returns diffusion map embedding for new data using the Nystrom extension
function pnew = diffMapRestrict(newData,evals,evecs,origData,eps)
dist = zeros(size(origData,1),1);
for iRow = 1:size(origData,1)
    dist(iRow) = norm(newData - origData(iRow,:));
end
w = exp(-(dist.^2)/eps^2);
k = (1/sum(w))*w;
pnew = (evecs' * k)./diag(evals);
end

