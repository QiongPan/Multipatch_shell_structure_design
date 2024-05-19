function [Ubar, Ctrlptsbar] = KnotInsertCurve(u_vec,U,Ctrlpts)

%%%%%% Insert knots for a NURBS curve %%%%%%
% input: 
%   u_vec      - knots to be inserted into U, a vector
%   U          - initial knot vector.(m+p+1)*1
%   Ctrlpts    - initial control points. (m+1)*4
% output: 
%   Ubar       - knot vector after insertion
%   Ctrlptsbar - control points after insertion
% Adapted from Algorithm A5.4 from 'The NURBS BOOK 2nd Edition' pg164.

[m,~] = size(Ctrlpts);
m = m - 1;
p = numel(U) - m -2;

u_vec = sort(u_vec);
U_ = [];

i = 1;
while(i <= numel(u_vec))
    u = u_vec(i);
    mult_ = FindMultiplicity(u, u_vec);
    i = i + mult_;
    mult = FindMultiplicity(u, U);
    if(mult_ > p - mult)
        mult_ = p - mult;
    end
    U_ = [U_,ones(1,mult_)*u];
end

u_vec = sort(U_);
r = numel(u_vec);
Ctrlptsbar = zeros(m+1+r, 4);

if (r ~= 0)    % insert r knots (Algorithm 5.4 insert r+1 knots)
    if(u_vec(1) < U(1) || u_vec(end) > U(end))
        error ('Trying to insert a knot outside the interval of definition')
    end
    
    Ubar = sort([U, u_vec]);
    a = FindSpan(m, p, u_vec(1), U);
    b = FindSpan(m, p, u_vec(end), U) + 1;
    Cp = zeros(m+1+r, 4);
    Cp(1:a-p+1, :) = Ctrlpts(1:a-p+1, :);
    Cp(b+r:end, :) = Ctrlpts(b:end, :);
    
    i = b + p - 1;  
    k = b + p + (r - 1);
    
    for j = (r-1) : -1 :0
        while(u_vec(j+1) <= U(i+1) && i > a)
            Cp(k-p, :) = Ctrlpts(i-p, :);
            k = k - 1;
            i = i - 1;
        end
        
        Cp(k-p, :) = Cp(k-p+1, :);
        
        for l = 1 : p
            ind = k - p + l;
            alpha = Ubar(k+l+1) - u_vec(j+1);
            if (abs(alpha) == 0.0)
                Cp(ind, :) = Cp(ind+1, :);
            else
                alpha = alpha / (Ubar(k+l+1) - U(i-p+l+1));
                Cp(ind, :) = alpha * Cp(ind, :) + (1.0-alpha) * Cp(ind+1, :);
            end          
        end
        
        k = k - 1;
    end
       
else 
    Cp = Ctrlpts;
    Ubar = U;    
end


end