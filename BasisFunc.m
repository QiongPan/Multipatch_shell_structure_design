/**
 * @file BasisFunc.m
 * @author Qiong Pan (PQ2019@mail.ustc.edu.cn)
 * @brief  compute the basis value for a B-Spline function.
 * Orginal work written by Qiong Pan in Matlab.
 * @date 2021-11-1
 */

function N = BasisFunc(i, p, u, U)

%%%%%% Basis function for B-Spline %%%%%%
% input: 
%   i          - the i-th knot span that u belongs to (obtained from FindSpan(),starting from 0) 
%   p          - the degree of spline functions
%   u          - parameter u \in [u_i, u_{i+1})
%   U          - knot vector  
% output: 
%   N          - N_j^k(u) j = i-p, ..., i. the values of spline basis functions at u, a vector
%              - dim = (p+1) * 1
% Adapted from Algorithm A2.2 from 'The NURBS BOOK 2nd Edition' pg70.

i =  i + 1;
N = zeros(p+1,1);
left = zeros(p+1,1);
right = zeros(p+1,1);

N(1) = 1.0;
for j = 1:p
    left(j+1) = u - U(i+1-j);
    right(j+1) = U(i+j) - u;
    saved = 0.0;

    for r=0:j-1
        temp = N(r+1)/(right(r+2) + left(j-r+1));
        N(r+1) = saved + right(r+2)*temp;
        saved = left(j-r+1)*temp;
    end

    N(j+1) = saved;
end

end
