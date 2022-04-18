function [ Eq, Dq ] = queue( eff_bw, Dq_max, lambda )

Dq = Dq_max;
a = 1/(eff_bw*Dq);
b = 1/(lambda*Dq);

x = 0:0.1:100;  % x = ln(1/Eq)
y1 = exp(a*x);
y2 = b*x + 1;

[xout,yout] = intersections(x,y1,x,y2,1);
x0 = max(xout);
Eq = exp(-x0);

end

