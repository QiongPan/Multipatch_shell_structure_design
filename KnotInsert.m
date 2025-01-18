/**
 * @file KnotInsert.m
 * @author Qiong Pan (PQ2019@mail.ustc.edu.cn)
 * @brief  knot insertion for a NURBS surface
 * Orginal work written by Qiong Pan in Matlab.
 * @date 2021-11-1
 */
function [Ubar, Vbar, Ctrlptsbar] = KnotInsert(u_vec,v_vec,U,V,Ctrlpts)

%%%%%% Insert u,v knots for a NURBS surface %%%%%%
% input: 
%   u_vec      - u knots to be inserted into UKnots, a vector
%   v_vec      - u knots to be inserted into UKnots, a vector
%   U          - initial u-knot vector.(m+p+1)*1
%   V          - initial v-knot vector.(n+q+1)*1
%   Ctrlpts    - initial control points. (m+1)*(n+1)*4(Homogeneous)
% output: 
%   Ubar       - u-knot vector after insertion
%   Vbar       - v-knot vector after insertion
%   Ctrlptsbar - control points after insertion
% Adapted from Algorithm A5.5 from 'The NURBS BOOK 2nd Edition' pg167.

[m,n,~] = size(Ctrlpts);
m = m - 1;  n = n - 1;
p = numel(U) - m -2;
q = numel(V) - n -2;

u_vec = sort(u_vec);
v_vec = sort(v_vec);
U_ = []; V_ = [];

i = 1;
while(i <= numel(u_vec))
    u = u_vec(i);
    mult_ = FindMultiplicity(u, u_vec);
    i = i + mult_;
    mult = FindMultiplicity(u, U);
    if(mult_ > p - mult)
        mult_ = p - mult;
    end
    U_ = [ones(1,mult_)*u, U_];
end

j = 1;
while(j <= numel(v_vec))
    v = v_vec(j);
    mult_ = FindMultiplicity(v, v_vec);
    j = j + mult_;
    mult = FindMultiplicity(v, V);
    if(mult_ > q - mult)
        mult_ = q - mult;
    end
    V_ = [ones(1,mult_)*v, V_];
end

u_vec = sort(U_);
v_vec = sort(V_);
r1 = numel(u_vec);
r2 = numel(v_vec);
Ctrlptsbar = zeros(m+1+r1, n+1+r2, 4);

if (r1 ~= 0)    
    if(u_vec(1) < U(1) || u_vec(end) > U(end))
        error ('Trying to insert a knot outside the interval of definition')
    end
    
    Ubar = sort([U, u_vec]);
    a = FindSpan(m, p, u_vec(1), U);
    b = FindSpan(m, p, u_vec(end), U) + 1;
    Cp = zeros(m+1+r1, n+1, 4);
    Cp(1:a-p+1, :, :) = Ctrlpts(1:a-p+1, :, :);
    Cp(b+r1:end, :, :) = Ctrlpts(b:end, :, :);
    
    i = b + p - 1;  
    k = b + p + r1 - 1;
    
    for j = r1-1 : -1 :0
        while(u_vec(j+1) <= U(i+1) && i > a)
            Cp(k-p, :, :) = Ctrlpts(i-p, :, :);
            k = k - 1;
            i = i - 1;
        end
        
        Cp(k-p, :, :) = Cp(k-p+1, :, :);
        
        for l = 1 : p
            ind = k - p + l;
            alpha = Ubar(k+l+1) - u_vec(j+1);
            if (abs(alpha) == 0.0)
                Cp(ind, :, :) = Cp(ind+1, :, :);
            else
                alpha = alpha / (Ubar(k+l+1) - U(i-p+l+1));
                Cp(ind, :, :) = alpha * Cp(ind, :, :) + (1.0-alpha) * Cp(ind+1, :, :);
            end          
        end
        
        k = k - 1;
    end
       
else 
    Cp = Ctrlpts;
    Ubar = U;    
end

if(r2 ~= 0)
    if(v_vec(1) < V(1) || v_vec(end) > V(end))
        error ('Trying to insert a knot outside the interval of definition')
    end
    
    Vbar = sort([V, v_vec]);
    a = FindSpan(n, q, v_vec(1), V);
    b = FindSpan(n, q, v_vec(end), V) + 1;    
    Ctrlptsbar = zeros(m+1+r1, n+1+r2, 4);
    Ctrlptsbar(:, 1:a-q+1, :) = Cp(:, 1:a-q+1, :);
    Ctrlptsbar(:, b+r2:end, :) = Cp(:, b:end, :);
    
    i = b + q - 1;  
    k = b + q + r2 - 1;
    
    for j = r2-1 : -1 : 0
        while(v_vec(j+1) <= V(i+1) && i > a)
            Ctrlptsbar(:, k-q, :) = Cp(:, i-q, :);
            k = k - 1;
            i = i - 1;
        end
        
        Ctrlptsbar(:, k-q, :) = Ctrlptsbar(:, k-q+1, :);
        
        for l = 1 : q
            ind = k - q + l;
            alpha = Vbar(k+l+1) - v_vec(j+1);
            if (abs(alpha) == 0.0)
                Ctrlptsbar(:, ind, :) = Ctrlptsbar(:, ind+1, :);
            else
                alpha = alpha / (Vbar(k+l+1) - V(i-q+l+1));
                Ctrlptsbar(:, ind, :) = alpha * Ctrlptsbar(:, ind, :) + (1.0-alpha) * Ctrlptsbar(:, ind+1, :);
            end          
        end
        
        k = k - 1;
        
    end
    
else   
    Ctrlptsbar = Cp;
    Vbar = V; 
end

end
