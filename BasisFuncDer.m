/**
 * @file BasisFuncDer.m
 * @author Qiong Pan (PQ2019@mail.ustc.edu.cn)
 * @brief compute 0-th to nth derivatives for Non-Zero B-Spline basis functions.
 * Orginal work written by Qiong Pan in Matlab.
 * @date 2021-11-1
 */

function NnthDer = BasisFuncDer(i, p, u, U, nth)

%%%%%% The 0-th to nth derivatives for Non-Zero B-Spline basis functions at u∈[u_i,u_{i+1}) %%%%%%
% input: 
%   i          - the i-th knot span that u belongs to (obtained from FindSpan(),starting from 0) 
%   p          - the degree of spline functions
%   u          - parameter u \in [u_i, u_{i+1})
%   U          - knot vector  
%   nth        - order of derivative(we get first to nth derivatives here) 
% output: 
%   NnthDer    - [N_r^p(u)]^{k} r = i-p, ..., i. k = 0,1,...nth, a matrix
%              - dimension = (nth+1) * (p+1)
% Adapted from Algorithm A2.3 from 'The NURBS BOOK 2nd Edition' pg72.

i = i+1; %% Matlab index of knots start from 1 not 0
ndu = zeros(p+1,p+1);
left = zeros(p+1,0);
right = zeros(p+1,0);
a = zeros(2,p+1);
ders = zeros(nth+1, p+1);

ndu(1,1) = 1.0;                  
for j = 1:p
	left(j+1) = u - U(i+1-j);
    right(j+1) = U(i+j) - u;
    saved = 0;
    for r = 0:j-1
        ndu(j+1,r+1) = right(r+2) + left(j-r+1);
        temp = ndu(r+1,j)/ndu(j+1,r+1);
        ndu(r+1,j+1) = saved + right(r+2)*temp;
        saved = left(j-r+1)*temp;
    end
    ndu(j+1,j+1) = saved;
end  

ders(1,1:p+1) = ndu(1:p+1,p+1);
    
for r = 0:p
    s1 = 0;
    s2 = 1;
    a(1,1) = 1;
    for k = 1:nth % compute kth derivative
        d = 0;
        rk = r-k;
        pk = p-k;
        if (r >= k)
            a(s2+1,1) = a(s1+1,1) / ndu(pk+2,rk+1);
            d = a(s2+1,1) * ndu(rk+1,pk+1);
        end
        if (rk >= -1)
            j1 = 1;
        else 
            j1 = -rk;
        end
        if (r-1 <= pk)
            j2 = k-1;
        else 
            j2 = p-r;
        end
        
        for j = j1:j2
            a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j))/ndu(pk+2,rk+j+1);
            d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1);
        end
        
        if (r <= pk)
            a(s2+1,k+1) = -a(s1+1,k) / ndu(pk+2,r+1);
            d = d + a(s2+1,k+1) * ndu(r+1,pk+1);
        end
        ders(k+1,r+1) = d;
        j = s1;
        s1 = s2;
        s2 = j;
    end    
end

r = p;
for k = 1:nth
    ders(k+1,1:p+1) = ders(k+1,1:p+1) * r;
    r = r * (p-k);
end
    
NnthDer = ders;

end
