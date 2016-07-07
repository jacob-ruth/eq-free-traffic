function [vec,val] = diffusionMap(epsilon,distMatrix,k)
    A = basicKernel(distMatrix);
    M = markovify(A);
    [eigenvec,eigenval] = eig(M);
    D = [1 0 0;
        0 3 0;
        0 0 2];
    P = [1 3 2;
        1 3 2;
        1 3 2]
    D2 = diag(sort(diag(D),'descend'))
    [c,ind] = sort(diag(D),'descend')
    P2 = P(:,ind)
    
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