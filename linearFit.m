function r_k = linearFit(evecs,k)
m = size(evecs,1);
PHI = evecs(:,1:k-1)'; %eigvectors are now rows
phi = evecs(:,k)';
ereg = median(pdist(PHI))/3;

minab = zeros(m,k);
for i = 1:m
    i
    minab(i,:) = fminsearch(@(params)toMin(params,phi,PHI,i),zeros(1,k));
end
alpha = minab(:,1);
beta = minab(:,2:end);

r_k2 = 0;
approx = zeros(1,m);
for i = 1:m
    approx(i) = alpha(i)+beta(i,:)*PHI(:,i);
    r_k2 = r_k2 + (phi(i)-(alpha(i)+beta(i,:)*PHI(:,i)))^2;
end
diff = phi - approx;
r_k = sqrt(r_k2/sum(phi.^2));


    function f = toMin(params,phi,P,i)
        a = params(1);
        b = params(2:end);
        f = 0;
        for j = 1:m
            if(j == i)
                continue;
            end
            dist = kernel(P(:,i),P(:,j));
            otherpart = (phi(j)- (a+b*P(:,j)))^2;
            
            f = f + dist*otherpart;
        end
    end

    function dist = kernel(r1,r2)
        dist = exp(-norm(r1-r2)^2/ereg^2);
    end
end