function [NewCtrlpts,NewU,NewV] = DegreeElevateSurface(p,q,U,V,Ctrlpts,tu,tv)

%%%%%% Degree elevate for a NURBS SURFACE %%%%%% 
% input:
%   p          - the degree of basis functions for u-curves.
%   q          - the degree of basis functions for v-curves.
%   Ctrlpts    - control points. 
%              - dim = (mc,nc,dim_c). 
%   U          - uknot vector
%   V          - vunot vector
%   tu         - Raise the Surface degree tu times in u-direction.
%   tv         - Raise the Surface degree tv times in v-direction.
% output:
%   NewCtrlpts - Control points after degree elevation
%   NewU       - uknot vector after degree elevation
%   NewV       - vknot vector after degree elevation
% Adapted from Algorithm A5.9 from 'The NURBS BOOK 2nd Edition' pg206.

[nc,mc,dim] = size(Ctrlpts);

%% u direction
if(tu > 0)        
    n = nc - 1;                                                           
    bezalfs = zeros(p+tu+1,p+1);     % coefficients for degree elevating the Bezier segments           
    bpts    = zeros(p+1,mc,dim);     % pth-degree Bezier control points of the current segment           
    ebpts   = zeros(p+tu+1,mc,dim);  % (p+tu)th-degree Bezier control points of the current segment           
    Nextbpts= zeros(p+1,mc,dim);     % leftmostpcontrol points of the next Bezier segment             
    alfs    = zeros(p,1);          	 % knot insertion alphas            

    m   = n + p + 1;                                  
    ph  = p + tu;                                           
    ph2 = floor(ph / 2);                                     

    %% compute bezier degree elevation coefficeients
    bezalfs(1,1) = 1;                                        
    bezalfs(ph+1,p+1) = 1;

    for i = 1:ph2                                           
        inv = 1/bincoeff(ph,i);                            
        mpi = min(p,i);                                   

        for j = max(0,i-tu):mpi                                 
            bezalfs(i+1,j+1) = inv*bincoeff(p,j)*bincoeff(tu,i-j); 
        end                                                       
    end                                                

    for i = ph2+1:ph-1                                    
        mpi = min(p,i);                                   
        for j = max(0,i-tu):mpi    
            bezalfs(i+1,j+1) = bezalfs(ph-i+1,p-j+1);          
        end                                                       
    end                                 

    mh = ph;                                  
    kind = ph + 1;                                       
    r = -1;                                                  
    a = p;                                  
    b = p + 1;                                 
    cind = 1;                                
    ua = U(1);                             
    NewCtrlpts(1,:,:) = Ctrlpts(1,:,:);
    NewU(1:ph+1) = ua;

    % initialize first bezier segment
    bpts(1:p+1,:,:) = Ctrlpts(1:p+1,:,:);
    
    %% big loop thru knot vector
    while (b < m)                                  
        i = b;                                                
        while (b < m && U(b+1) == U(b+2))                       
            b = b + 1;                                         
        end                                                   
        mul = b - i + 1;                                      
        mh = mh + mul + tu;                                     
        ub = U(b+1);                                        
        oldr = r;                                              
        r = p - mul;                                          

       % insert knot u(b) r times
        if oldr > 0                                         
            lbz = floor((oldr+2)/2);                          
        else                                                
            lbz = 1;                                            
        end                                                   

        if r > 0                                            
            rbz = ph - floor((r+1)/2);                          
        else                                                  
            rbz = ph;                                           
        end

        % insert knot to get bezier segment
        if r > 0
            numer = ub - ua;                                  
            for k = p:-1:mul+1
                alfs(k-mul) = numer / (U(a+k+1)-ua);            
            end                                           

            for j = 1:r                                        
                save = r - j;                               
                s = mul + j;                              

                for k = p:-1:s      
                    bpts(k+1,:,:) = alfs(k-s+1) * bpts(k+1,:,:) + (1-alfs(k-s+1)) *  bpts(k,:,:);                                            
                end   

                Nextbpts(save+1,:,:) = bpts(p+1,:,:);                                               
            end                                             
        end   % end of insert knot

        %% degree elevate bezier(only points lbz,...,ph are used below)
        for i = lbz:ph 
            ebpts(i+1,:,:) = 0;                                                  
            mpi = min(p, i); 

            for j = max(0,i-tu):mpi    
                ebpts(i+1,:,:) = ebpts(i+1,:,:) + bezalfs(i+1,j+1) * bpts(j+1,:,:);                                                
            end                                                    
        end  % end of degree elevating bezier

        %% must remove knot u = U[a] oldr times
        if (oldr > 1)   
            first = kind - 2;                    
            last = kind;                         
            den = ub - ua;                       
            bet = floor((ub-NewU(kind)) / den);     

            % knot removal loop
            for tr = 1:oldr-1                                    
                i = first;                                      
                j = last;                                        
                kj = j - kind + 1;

                while j-i > tr   % loop and compute the new control points
                if i < cind                                  
                    alf = (ub-NewU(i+1)) / (ua-NewU(i+1));
                    NewCtrlpts(i+1,:,:) = alf*NewCtrlpts(i+1,:,:) + (1-alf)*NewCtrlpts(i,:,:);                                         
                end                                       
                if j >= lbz                                 
                    if j-tr <= kind-ph+oldr                    
                        gam = (ub - NewU(j-tr+1)) / den; 
                        ebpts(kj+1,:,:) = gam*ebpts(kj+1,:,:) + (1-gam) * ebpts(kj+2,:,:);                                                                                  
                    else     
                        ebpts(kj+1,:,:) = bet*ebpts(kj+1,:,:) + (1-bet) * ebpts(kj+2,:,:);                                      
                    end                                      
                end                                         
                i = i + 1;                                
                j = j - 1;                                 
                kj = kj - 1;                             
                end                                         

                first = first - 1;                             
                last = last + 1;                             
            end                                                 
        end   % end of removing knot u=U[a]

        %%
        if a ~= p           % load the knot ua                                   
            for i = 0:ph-oldr-1                                   
                NewU(kind+1) = ua;                                 
                kind = kind + 1;                                
            end
        end                                                   

        for j = lbz:rbz    	% load ctrlpts into NewCtrlpts                             
            NewCtrlpts(cind+1,:,:) = ebpts(j+1,:,:);                                               
            cind = cind + 1;                           
        end                                         

        if b < m            % setup for next pass thru loop
            bpts(1:r,:,:) = Nextbpts(1:r,:,:);
            bpts(r+1:p+1,:,:) = Ctrlpts(b-p+r+1:b+1,:,:);                                               
            a = b;                                          
            b = b + 1;                                       
            ua = ub;                                         
        else   
            NewU(kind+(0:ph)+1) = ub;                                              
        end  

    end   % end of while-loop(b < m)

else
    NewCtrlpts = Ctrlpts;
    NewU = U;
end

%% v direction
Ctrlpts = NewCtrlpts;
[nc,mc,dim] = size(Ctrlpts);
if (tv > 0)
    n = mc - 1;                                                           
    bezalfs = zeros(q+tv+1,q+1);   	% coefficients for degree elevating the Bezier segments           
    bpts    = zeros(nc,q+1,dim);    % qth-degree Bezier control points of the current segment           
    ebpts   = zeros(nc,q+tv+1,dim); % (q+tv)th-degree Bezier control points of the current segment           
    Nextbpts= zeros(nc,q+1,dim);    % leftmostpcontrol points of the next Bezier segment             
    alfs    = zeros(q,1);          	% knot insertion alphas            

    m   = n + q + 1;                                  
    ph  = q + tv;                                           
    ph2 = floor(ph / 2);                                     

    % compute bezier degree elevation coefficeients
    bezalfs(1,1) = 1;                                        
    bezalfs(ph+1,q+1) = 1;

    for i = 1:ph2                                           
        inv = 1/bincoeff(ph,i);                            
        mpi = min(q,i);                                   

        for j = max(0,i-tv):mpi                                 
            bezalfs(i+1,j+1) = inv*bincoeff(q,j)*bincoeff(tv,i-j); 
        end                                                       
    end                                                

    for i = ph2+1:ph-1                                    
        mpi = min(q,i);                                   
        for j = max(0,i-tv):mpi    
            bezalfs(i+1,j+1) = bezalfs(ph-i+1,q-j+1);          
        end                                                       
    end                                 

    mh = ph;                                  
    kind = ph + 1;                                       
    r = -1;                                                  
    a = q;                                  
    b = q + 1;                                 
    cind = 1;                                
    va = V(1);                             
    NewCtrlpts(:,1,:) = Ctrlpts(:,1,:);
    NewV(1:ph+1) = va;

    % initialise first bezier seg
    bpts(:,1:q+1,:) = Ctrlpts(:,1:q+1,:);
    while (b < m)  % big loop thrv knot vector                                
        i = b;                                                
        while (b < m && V(b+1) == V(b+2))                       
            b = b + 1;                                         
        end                                                   
        mul = b - i + 1;                                      
        mh = mh + mul + tv;                                     
        vb = V(b+1);                                        
        oldr = r;                                              
        r = q - mul;                                          

       % insert knot v(b) r times
        if oldr > 0                                         
            lbz = floor((oldr+2)/2);                          
        else                                                
            lbz = 1;                                            
        end                                                   

        if r > 0                                            
            rbz = ph - floor((r+1)/2);                          
        else                                                  
            rbz = ph;                                           
        end

        % insert knot to get bezier segment
        if r > 0
            numer = vb - va;                                  
            for k = q:-1:mul+1
                alfs(k-mul) = numer / (V(a+k+1)-va);            
            end                                           

            for j = 1:r                                        
                save = r - j;                               
                s = mul + j;                              

                for k = q:-1:s      
                    bpts(:,k+1,:) = alfs(k-s+1) * bpts(:,k+1,:) + (1-alfs(k-s+1)) *  bpts(:,k,:);                                            
                end   

                Nextbpts(:,save+1,:) = bpts(:,q+1,:);                                               
            end                                             
        end   % end of insert knot

        % degree elevate bezier(only points lbz,...,ph are used below)
        for i = lbz:ph 
            ebpts(:,i+1,:) = 0;                                                  
            mpi = min(q, i); 

            for j = max(0,i-tv):mpi    
                ebpts(:,i+1,:) = ebpts(:,i+1,:) + bezalfs(i+1,j+1) * bpts(:,j+1,:);                                                
            end                                                    
        end  % end of degree elevating bezier

        % must remove knot v = V[a] oldr times
        if (oldr > 1)   
            first = kind - 2;                    
            last = kind;                         
            den = vb - va;                       
            bet = floor((vb-NewV(kind)) / den);     

            % knot removal loop
            for tr = 1:oldr-1                                    
                i = first;                                      
                j = last;                                        
                kj = j - kind + 1;

                while j-i > tr   % loop and compute the new control points
                if i < cind                                  
                    alf = (vb-NewV(i+1)) / (va-NewV(i+1));
                    NewCtrlpts(:,i+1,:) = alf*NewCtrlpts(:,i+1,:) + (1-alf)*NewCtrlpts(:,i,:);                                         
                end                                       
                if j >= lbz                                 
                    if j-tr <= kind-ph+oldr                    
                        gam = (vb - NewV(j-tr+1)) / den; 
                        ebpts(:,kj+1,:) = gam*ebpts(:,kj+1,:) + (1-gam) * ebpts(:,kj+2,:);                                                                                  
                    else     
                        ebpts(:,kj+1,:) = bet*ebpts(:,kj+1,:) + (1-bet) * ebpts(:,kj+2,:);                                      
                    end                                      
                end                                         
                i = i + 1;                                
                j = j - 1;                                 
                kj = kj - 1;                             
                end                                         

                first = first - 1;                             
                last = last + 1;                             
            end                                                 
        end   % end of removing knot v=V[a]

        if a ~= q           % load the knot va                                   
            for i = 0:ph-oldr-1                                   
                NewV(kind+1) = va;                                 
                kind = kind + 1;                                
            end
        end                                                   

        for j = lbz:rbz    	% load ctrlpts into NewCtrlpts                             
            NewCtrlpts(:,cind+1,:) = ebpts(:,j+1,:);                                               
            cind = cind + 1;                           
        end                                         

        if b < m            % setup for next pass thrv loop
            bpts(:,1:r,:) = Nextbpts(:,1:r,:);
            bpts(:,r+1:q+1,:) = Ctrlpts(:,b-q+r+1:b+1,:);                                               
            a = b;                                          
            b = b + 1;                                       
            va = vb;                                         
        else   
            NewV(kind+(0:ph)+1) = vb;                                              
        end  

    end   % end of while-loop(b < m)

else
    NewCtrlpts = Ctrlpts;
    NewV = V;   
end

end


function b = bincoeff(n,k)
%  Computes the binomial coefficient.
%
%      ( n )      n!
%      (   ) = --------
%      ( k )   k!(n-k)!
%  b = bincoeff(n,k)
%
%  Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.
                                                        
b = floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));    
end                                                      


function f = factln(n)
% computes ln(n!)
if n <= 1, f = 0; return, end
f = gammaln(n+1); %log(factorial(n));
end


%%
% function [NewCtrlpts,NewU] = DegreeElevateCurve(p,U,Ctrlpts,t)
% 
% % Degree elevation for a NURBS CURVE. 
% % input:
% %   p           - Degree of the NURBS curve.
% %   Ctrlpts     - Control points. dim = (nc,dim).
% %   U           - Knot vector
% %   t           - Raise the B-Spline degree t times.
% % output:
% %   NewCtrlpts  - Control points of the new curve. 
% %   NewU        - Knot vector of the new curve. 
% % Adapted from Algorithm A5.9 from 'The NURBS BOOK 2nd Edition' pg70.
% 
% 
% [nc,dim] = size(Ctrlpts);
% n = nc - 1;                                       
%                                  
% bezalfs =  zeros(p+t+1,p+1);    % coefficients for degree elevating the Bezier segments           
% bpts    = zeros(p+1,dim);      	% pth-degree Bezier control points of the current segment           
% ebpts   = zeros(p+t+1,dim);     % (p+t)th-degree Bezier control points of the current segment           
% Nextbpts= zeros(p+1,dim);      	% leftmostpcontrol points of the next Bezier segment             
% alfs    = zeros(p,1);          	% knot insertion alphas            
%    
% m   = n + p + 1;                                  
% ph  = p + t;                                           
% ph2 = floor(ph / 2);                                     
%                                                         
% % compute bezier degree elevation coefficeients
% bezalfs(1,1) = 1;                                        
% bezalfs(ph+1,p+1) = 1;
% 
% for i = 1:ph2                                           
% 	inv = 1/bincoeff(ph,i);                            
% 	mpi = min(p,i);                                   
%                                                        
%     for j = max(0,i-t):mpi                                 
%         bezalfs(i+1,j+1) = inv*bincoeff(p,j)*bincoeff(t,i-j); 
%     end                                                       
% end                                                
%                                                    
% for i = ph2+1:ph-1                                    
% 	mpi = min(p,i);                                   
%     for j = max(0,i-t):mpi    
%         bezalfs(i+1,j+1) = bezalfs(ph-i+1,p-j+1);          
%     end                                                       
% end                                 
%                                      
% mh = ph;                                  
% kind = ph + 1;                                       
% r = -1;                                                  
% a = p;                                  
% b = p + 1;                                 
% cind = 1;                                
% ua = U(1);                             
% NewCtrlpts(1,:) = Ctrlpts(1,:);
% NewU(1:ph+1) = ua;
%  
% % initialise first bezier seg
% bpts(1:p+1,:) = Ctrlpts(1:p+1,:);
% while (b < m)  % big loop thru knot vector                                
%     i = b;                                                
%     while (b < m && U(b+1) == U(b+2))                       
%         b = b + 1;                                         
%     end                                                   
%     mul = b - i + 1;                                      
%     mh = mh + mul + t;                                     
%     ub = U(b+1);                                        
%     oldr = r;                                              
%     r = p - mul;                                          
%                                                           
%    % insert knot u(b) r times
%     if oldr > 0                                         
%         lbz = floor((oldr+2)/2);                          
%     else                                                
%         lbz = 1;                                            
%     end                                                   
%    
% 	if r > 0                                            
%         rbz = ph - floor((r+1)/2);                          
%     else                                                  
%         rbz = ph;                                           
%     end
%     
%     % insert knot to get bezier segment
%     if r > 0
%         numer = ub - ua;                                  
%         for k = p:-1:mul+1
%             alfs(k-mul) = numer / (U(a+k+1)-ua);            
%         end                                           
%       
%         for j = 1:r                                        
%             save = r - j;                               
%             s = mul + j;                              
%                                                       
%             for k = p:-1:s      
%                 bpts(k+1,:) = alfs(k-s+1) * bpts(k+1,:) + (1-alfs(k-s+1)) *  bpts(k,:);                                            
%             end   
%          
%             Nextbpts(save+1,:) = bpts(p+1,:);                                               
%         end                                             
% 	end   % end of insert knot
% 
%     % degree elevate bezier(only points lbz,...,ph are used below)
% 	for i = lbz:ph 
%         ebpts(i+1,:) = 0;                                                  
%         mpi = min(p, i); 
%       
%         for j = max(0,i-t):mpi    
%             ebpts(i+1,:) = ebpts(i+1,:) + bezalfs(i+1,j+1) * bpts(j+1,:);                                                
%         end                                                    
%     end  % end of degree elevating bezier
%     
%     % must remove knot u = U[a] oldr times
% 	if oldr > 1   
%         first = kind - 2;                    
%         last = kind;                         
%         den = ub - ua;                       
%         bet = floor((ub-NewU(kind)) / den);     
%         
%         % knot removal loop
%         for tr = 1:oldr-1                                    
%             i = first;                                      
%             j = last;                                        
%             kj = j - kind + 1;
%          
%             while j-i > tr   % loop and compute the new control points
%             if i < cind                                  
%                 alf = (ub-NewU(i+1)) / (ua-NewU(i+1));
%                 NewCtrlpts(i+1,:) = alf*NewCtrlpts(i+1,:) + (1-alf)*NewCtrlpts(i,:);                                         
%             end                                       
%          	if j >= lbz                                 
%              	if j-tr <= kind-ph+oldr                    
%                     gam = (ub - NewU(j-tr+1)) / den; 
%                     ebpts(kj+1,:) = gam*ebpts(kj+1,:) + (1-gam) * ebpts(kj+2,:);                                                                                  
%                 else     
%                     ebpts(kj+1,:) = bet*ebpts(kj+1,:) + (1-bet) * ebpts(kj+2,:);                                      
%                 end                                      
%             end                                         
%             i = i + 1;                                
%             j = j - 1;                                 
%             kj = kj - 1;                             
%             end                                         
%                                                          
%           	first = first - 1;                             
%             last = last + 1;                             
%         end                                                 
%     end   % end of removing knot u=U[a]
% 	
% 	if a ~= p           % load the knot ua                                   
%         for i = 0:ph-oldr-1                                   
%             NewU(kind+1) = ua;                                 
%             kind = kind + 1;                                
%         end
%     end                                                   
%    
%     for j = lbz:rbz    	% load ctrlpts into NewCtrlpts                             
%     	NewCtrlpts(cind,:) = ebpts(j+1,:);                                               
%        	cind = cind + 1;                           
%     end                                         
%                                   
%   	if b < m            % setup for next pass thru loop
%         bpts(1:r,:) = Nextbpts(1:r,:);
%         bpts(r+1:p+1,:) = Ctrlpts(b-p+r+1:b+1,:);                                               
%         a = b;                                          
%         b = b + 1;                                       
%         ua = ub;                                         
%     else   
%         NewU(kind+(0:ph)+1) = ub;                                              
%     end  
% 
% end   % end of while-loop(b < m)
% 
% end  
                                                         