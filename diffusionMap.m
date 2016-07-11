function [vec,val] = diffusionMap(epsilon,distMatrix,k)
    %find Markov matrix and its eigenvalues/eigenvectors
    A = basicKernel(distMatrix);
    M = markovify(A);
    [eigenvec,eigenval] = eig(M);
    
    %sort eigenvalues and eigenvectors
    seval = diag(sort(diag(eigenval),'descend'));
    [c,ind] = sort(diag(eigenval),'descend');
    sevec = eigenvec(:,ind);
    seval = sparse(seval);
    
    %return first k eigenvalues/vectors
    val = seval(2:k+1,2:k+1);
    vec = sevec(:,2:k+1);
    
    function M = markovify(af)
        M = zeros(size(af));
        rsum = sum(af,2);
        for r = 1:length(rsum)
            M(r,:) = af(r,:)./rsum(r);
        end
    end

    function dist = basicKernel(s)
        dist = exp(-s/epsilon^2);
    end
end