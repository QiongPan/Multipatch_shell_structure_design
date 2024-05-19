function r = FindSpan(n, p, u, U)

% Determine the knot span index             
% input:
%   n          - number of control points - 1, starting from 0
%   p          - spline degree
%   u          - parameter
%   U          - knot vector
% output:
%   r          - knot span index, u \in [u_r, u_{r+1})  p<=r<=n

DefDomain_U = U(p+1:n+2);
if (u > DefDomain_U(end) || u < DefDomain_U(1))
    error('u is out of range.')
end

% if u == U(end)
%     r = n; 
%     return;
% end

r = find(u >= DefDomain_U,1,'last') + p - 1;    % find(X,n,direction),finds the last n indices
                                                % corresponding to nonzero elements in X. 

if r > n
    r = n; 
end

end