function f=frosenbrock(x)
    %if size(x,1) < 2 error('dimension must be greater one'); end
    %f = 100*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);
    f = norm(x);
    %f = (1.5 - x(1) + x(1)*x(2)).^2 + (2.25 - x(1) + x(1)*(x(2).^2)).^2 ...
        + (2.625 - x(1) + x(1)*(x(2).^3)).^2;