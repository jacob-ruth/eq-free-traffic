function checkNewton()
    w = [3;1];
    newGuess = [.25;-.75];
    ref_2 = 0;
    
    dx = .001;
    dy = .001;

    u = newGuess;
    for k=1:10
        fprintf('starting iteration %f \n', k);
        f = F(ref_2, u(2),u(1));
        fprintf('f is %f \n', f)
        Df = jacobian(ref_2, u(1), u(2), w);
        neww = w(1)*(u(1)-newGuess(1)) + w(2)*(u(2) - newGuess(2))
        u = u - Df^(-1)*[f;w(1)*(u(1)-newGuess(1)) + w(2)*(u(2) - newGuess(2))]
    end

    function f = F(~, y, x)
        f = 1-(x^2 + 2*y^2);
    end

    function j = jacobian(~, x, y, w)
%         j = [-2*x, -4*y;w'];
        j(1,1) = (F(0,y,x+dx)-F(0,y,x))/dx;
        j(1,2) = (F(0,y+dy,x)-F(0,y,x))/dy;
        j(2,:) = w';
    end
end