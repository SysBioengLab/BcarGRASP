function [fc,grad] = detailedEnzymeConstrained(x,var1,var2,sign)

fc = sign.*x(var1)-sign.*x(var2);
if (nargout > 1)
    grad = zeros(1,51);
    grad(1,var1) = sign;
    grad(1,var2) = -sign;
end
end