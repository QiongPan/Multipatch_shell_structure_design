function [s, w] = GaussInt(a,b,ngp)

%%%%%% Guass points and weights for integration %%%%%%
% input: 
%   a          - left limits of integration span
%   b          - right limits of integration span
%   ngp        - number of guass points
% output: 
%   s          - guass point vector
%   w          - weight vector

switch(ngp)
    case 1
        s = 0;
        w = 2;
    case 2
        s = [-sqrt(1/3), sqrt(1/3)]';
        w = [1, 1]';
    case 3
        s = [-0.774596669241, 0, 0.774596669241]';
        w = [5/9, 8/9, 5/9]';
    case 4
        s = [-0.861136311594053, -0.339981043584856, 0.339981043584856,0.861136311594053]';
        w = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]';
    case 5
        s = [-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664]';
        w = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189]';
end

% t \in [-1,1], x = t*(b-a)/2+(a+b)/2 \in [a,b]
h = (b-a)/2;
m = (a+b)/2;
s = s * h + m;
w = w * h;  % coefficients when converting integrals: dx/dt

end