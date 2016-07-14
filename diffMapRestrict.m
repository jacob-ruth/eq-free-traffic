%% restriction operator based on diffusion map
% newData   - point to get restriction values for (column vector)
% evals     - square matrix of eigenvalues
% evecs     - matrix with eigenvectors in columns
% origData  - matrix of data in columns used to compute eigenvalues/vectors
% eps       - length scale for kernel
% 
% returns diffusion map embedding for new data using the Nystrom extension
% in a column vector
function pnew = diffMapRestrict(newData,evals,evecs,origData,eps)
dist = zeros(size(origData,2),1);
for iCol = 1:size(origData,2)
    dist(iCol) = norm(newData - origData(:,iCol));
end
w = exp(-(dist.^2)/eps^2);
k = (1/sum(w))*w;
pnew = (evecs' * k)./diag(evals);
end

