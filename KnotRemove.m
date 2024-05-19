function [Ubar, Vbar, Ctrlptsbar, tu, tv] = KnotRemove(U,V,Ctrlpts,u,v,num_u,num_v,TOL)

%%%%%% Remove knot u num_u times and knot v num_v times for a NURBS surface %%%%%%
% input:
%   U          - initial u-knot vector UKnots.(m+p+1)*1
%   V          - initial v-knot vector VKnots.(n+q+1)*1
%   Ctrlpts    - initial control points (m+1) * (n+1) * 4
%   u          - u knot to be removed from U.
%   v          - v knot to be removed from V.
%   num_u      - remove u from U for num_u times.
%   num_v      - remove v from V for num_v times.
%   TOL        - tolerance for deviation before and after removal
% output:
%   Ubar       - u-knots after removal.
%   Vbar       - v-knots after removal.
%   Ctrlptsbar - control points after removal.
%   tu         - actual times for u been removed from U
%   tv         - actual times for v been removed from V 
% Adapted from Algorithm A5.8 from 'The NURBS BOOK 2nd Edition' pg185

[m,n,~] = size(Ctrlpts);
m = m - 1;  n = n - 1;
p = numel(U) - m -2;
q = numel(V) - n -2;
ru = FindSpan(m, p, u, U);
rv = FindSpan(n, q, v, V);
su = FindMultiplicity(u, U);
sv = FindMultiplicity(v, V);

[Ubar, CtrlptsUbar, tu] = UKnotRemove(m,p,U,Ctrlpts,u,ru,su,num_u,TOL);
[Vbar, Ctrlptsbar, tv] = VKnotRemove(n,q,V,CtrlptsUbar,v,rv,sv,num_v,TOL);

end

function [Ubar, Ctrlptsbar, t] = UKnotRemove(m,p,U,Pw,u,r,s,num, TOL)

% Remove knot u num times to obtain new control points and knots.
% input: m,p,U,Ctrlpts,u,r,s,num
%       m  - maximum index of control points in u-direction. (count = m+1)
%       p  - degree in u direction.
%       U  - initial u knots vector UKnots.(m+p+1)*1.
%       Ctrlpts - initial control points (m+1)*(n+1)*4.
%       u  - u knot to be removed from U.
%       r  - u = u_r ≠ u_{r+1}. u-span where u lied.
%       num- remove u num times.
% output: Ubar, Ctrlptsbar, t 
%       Ubar  - new u-knots after removal.
%       Ctrlptsbar  - new control points after removal.
%       t     - actual times for u been removed from U. 
% Adapted from Algorithm A5.5 from 'The NURBS BOOK 2nd Edition' pg167.

if(num > s)
    num = s;
end

ord = p + 1;
fout = floor((2*r-s-p)/2); % First control point out
last = r - s;
first = r - p;
temp = zeros(last-first+1,size(Pw,2), 4);

t = 0;
while(t < num)
    off = first - 1;    % diff in index between temp and P
    temp(1,:,:) = Pw(off+1,:,:); 
    temp(last+2-off,:,:) = Pw(last+2,:,:);
    i = first; 
    j = last;
    ii = 1; 
    jj = last - off;
    remflag = 0;
    
    % Compute new control points for one removal step
    while(j - i > t)        
        alfi = (u - U(i+1)) / (U(i+ord+t+1) - U(i+1));
        alfj = (u - U(j-t+1)) / (U(j+ord+1) - U(j-t+1));
        temp(ii+1,:,:) = (Pw(i+1,:,:)-(1-alfi)*temp(ii,:,:)) / alfi;
        temp(jj+1,:,:) = (Pw(j+1,:,:)-alfj*temp(jj+2,:,:)) / (1-alfi);
        i = i + 1;  ii = ii + 1;
        j = j - 1;  jj = jj - 1;
    end
    
    % check if knot removable
    if (j - i < t)  
        if(Distance4D(temp(ii,:,:),temp(jj+2,:,:)) <= TOL)
            remflag = 1;
        end        
    else
        alfi = (u - U(i+1)) / (U(i+ord+t+1)-U(i+1));
        if(Distance4D(Pw(i+1,:,:),alfi*temp(ii+t+2,:,:)+(1-alfi)*temp(ii,:,:)) <= TOL)
            remflag = 1;
        end
    end  
    
    if (remflag == 0)   % Cannot remove any more knots
        
        break;          % get out of for-loop
    else                % Sunccessful removal. Save new control points
        i = first;  j = last;
        
        while (j - i > t)
            Pw(i+1,:,:) = temp(i-off+1,:,:);
            Pw(j+1,:,:) = temp(j-off+1,:,:);
            i = i + 1;  j = j - 1;            
        end
         
    end
    
    first = first - 1;
    last = last + 1;   
    
    t = t + 1;
end

if (t == 0) 
    Ubar = U(1:end-t);
    Ctrlptsbar = Pw(1:end-t,:,:);
    return;
end

% shift knots;
for k = r+1:m+p+1
    U(k-t+1) = U(k+1);   
end
Ubar = U(1:end-t);

% Pj thru Pi will be overwritten
j = fout;
i = j;
for k = 1:t-1
    if(mod(k,2) == 1)
        i = i + 1;
    else
        j = j - 1;
    end
end

for k = i+1:m
    Pw(j+1,:,:) = Pw(k+1,:,:);
    j = j + 1;
end

Ctrlptsbar = Pw(1:end-t,:,:);

end

function [Vbar, Ctrlptsbar, t] = VKnotRemove(n,q,V,Pw,v,r,s,num, TOL)

if(num > s)
    num = s;
end

ord = q + 1;
fout = floor((2*r-s-q)/2); % First control point out
last = r - s;
first = r - q;
temp = zeros(size(Pw,1), last-first+1,4);

t = 0;
while(t < num)
    off = first - 1;    % diff in index between temp and P
    temp(:,1,:) = Pw(:,off+1,:); 
    temp(:,last+2-off,:) = Pw(:,last+2,:);
    i = first; 
    j = last;
    ii = 1; 
    jj = last - off;
    remflag = 0;
    
    % Compute new control points for one removal step
    while(j - i > t)        
        alfi = (v - V(i+1)) / (V(i+ord+t+1) - V(i+1));
        alfj = (v - V(j-t+1)) / (V(j+ord+1) - V(j-t+1));
        temp(:,ii+1,:) = (Pw(:,i+1,:)-(1-alfi)*temp(:,ii,:)) / alfi;
        temp(:,jj+1,:) = (Pw(:,j+1,:)-alfj*temp(:,jj+2,:)) / (1-alfi);
        i = i + 1;  ii = ii + 1;
        j = j - 1;  jj = jj - 1;
    end
    
    % check if knot removable
    if (j - i < t)  
        if(Distance4D(temp(:,ii,:),temp(:,jj+2,:)) <= TOL)
            remflag = 1;
        end        
    else
        alfi = (v - V(i+1)) / (V(i+ord+t+1)-V(i+1));
        if(Distance4D(Pw(:,i+1,:),alfi*temp(:,ii+t+2,:)+(1-alfi)*temp(:,ii,:)) <= TOL)
            remflag = 1;
        end
    end  
    
    if (remflag == 0)   % Cannot remove any more knots
        break;          % get out of for-loop
    else                % Sunccessful removal. Save new control points
        i = first;  j = last;
        
        while (j - i > t)
            Pw(:,i+1,:) = temp(:,i-off+1,:);
            Pw(:,j+1,:) = temp(:,j-off+1,:);
            i = i + 1;  j = j - 1;            
        end
         
    end
    
    first = first - 1;
    last = last + 1;   
    
    t = t + 1;
end

if (t == 0) 
    Vbar = V(1:end-t);
    Ctrlptsbar = Pw(1:end-t,:,:);
    return;
end

% shift knots;
for k = r+1:n+q+1
    V(k-t+1) = V(k+1);   
end
Vbar = V(1:end-t);

% Pj thru Pi will be overwritten
j = fout;
i = j;
for k = 1:t-1
    if(mod(k,2) == 1)
        i = i + 1;
    else
        j = j - 1;
    end
end

for k = i+1:n
    Pw(:,j+1,:) = Pw(:,k+1,:);
    j = j + 1;
end

Ctrlptsbar = Pw(:,1:end-t,:);

end

function dist4D = Distance4D(Pw, Qw)

    dist4D = 0;
    [dim1,dim2,~] = size(Pw);
    for i = 1:dim1
        for j = 1:dim2
            P = reshape(Pw(i,j,:),1,[]);
            Q = reshape(Qw(i,j,:),1,[]);
            dist4D = max(norm(P-Q),dist4D);
        end
    end
    
end

%% rational 和 non-rational 的deviation调整和TOL选择