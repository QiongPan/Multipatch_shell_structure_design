function F = ComputeF(mAnalysis,nAnalysis,p,q,UAnalysis,VAnalysis,UEle,VEle,Ctrlpts,ploadsMat,loadType, ne_u, ne_v)

%%%%%% Compute the force vector in correlance to control points, the force is applied on mid-surface and integrated along thickness direction %%%%%%
% input:
%   mAnalysis  - number of control points - 1 in analysis model for u-curves (U-direction)/(s-direction)
%   nAnalysis  - number of control points - 1 in analysis model for v-curves (V-direction)/(t-direction)
%   p          - spline degree                          (U-direction)
%   q          - spline degree                          (V-direction)
%   UAnalysis  - knot vector in analysis model          (U-direction)
%   VAnalysis  - knot vector in analysis model          (V-direction)
%   UEle       - unique uknots for knot spans           (U-direction)
%   VEle       - unique vknots for knot spans           (V-direction)
%   Ctrlpts    - control points
%   ploads_Mat - the load applied on the patch
%              - a matrix with 5 columns
%              - ploads_Mat(i,:) = (s_i,t_i,fx_i,fy_i,fz_i)
%   loadType   - 0:point force; 
%              - 1-3: line force
%   ne_u       - number of element(knot spans) in u-direction    (U-direction)
%   ne_v       - number of element(knot spans) in v-direction    (V-direction)
% output:
%   F          - global force vector (K*U = F)

global ng_u ng_v ng_th thickness;
F = zeros(5*(mAnalysis+1)*(nAnalysis+1), 1);      
[s3_vec, w3_vec] = GaussInt(-1,1,ng_th);

switch loadType
    case 0 % point force Fe = Me'* P  P = [fx,fy,fz] at point (s,t)
        Weights = reshape(Ctrlpts(:,:,4), mAnalysis+1,nAnalysis+1);
        
        for kk = 1:size(ploadsMat,1)
            s0 = ploadsMat(kk,1);
            t0 = ploadsMat(kk,2);
            pload = ploadsMat(kk,3:5);
            
            dd = 0.1/5;
            uu = (s0-0.05):dd:(s0+0.05); uu = uu((uu>=0)); uu = uu((uu<=1));
            vv = (t0-0.05):dd:(t0+0.05); vv = vv((vv>=0)); vv = vv((vv<=1));
            
            % just for s5
            uu = [s0];
            vv = [t0];
            % just for s5
            
            
            pload = pload / numel(uu) /numel(vv);
            for i1 = 1:numel(uu)
                s0 = uu(i1);
                for j1 = 1:numel(vv)
                    t0 = vv(j1);
                    [IdxU, IdxV, N, Ns, Nt, ~, ~, ~] = NrbSurDer2(mAnalysis,nAnalysis,p,q,s0,t0,UAnalysis,VAnalysis,Weights);
                    PIJ = Ctrlpts(IdxU+1,IdxV+1,1:3);
                    Ss = reshape(sum(PIJ.*Ns, [1,2]),1,3);
                    St = reshape(sum(PIJ.*Nt, [1,2]),1,3);
                    
                    % local coordinate axises
                    v3 = cross(Ss,St);
                    v3 = v3 / norm(v3);
                    St_norm = norm(St);
                    Ss_norm = norm(Ss);
                    if(Ss_norm)
                        v1 = Ss/Ss_norm;
                        v2 = cross(v3, v1);
                    elseif(St_norm)
                        v2 = St/St_norm;
                        v1 = cross(v2, v3);
                    end
                    l1 = v1(1); m1 = v1(2); n1 = v1(3);
                    l2 = v2(1); m2 = v2(2); n2 = v2(3);
                    
                    NIJ  = reshape(N,  1, []);
                    [J, I] = meshgrid(IdxV, IdxU);
                    Idx = sub2ind([mAnalysis+1, nAnalysis+1], I(:)+1, J(:)+1);
                    dispIdxIJ = reshape(((Idx-1)*5+(1:5))',1,[]);
                    
                    
                    for a = 1:numel(w3_vec)
                        theta = s3_vec(a);
                        M = zeros(3,5*(p+1)*(q+1));
                        M(1,1:5:5*(p+1)*(q+1)) =  NIJ;
                        M(1,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*l2*NIJ;
                        M(1,5:5:5*(p+1)*(q+1)) =  theta*thickness/2*l1*NIJ;
                        M(2,2:5:5*(p+1)*(q+1)) =  NIJ;
                        M(2,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*m2*NIJ;
                        M(2,5:5:5*(p+1)*(q+1)) =  theta*thickness/2*m1*NIJ;
                        M(3,3:5:5*(p+1)*(q+1)) =  NIJ;
                        M(3,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*n2*NIJ;
                        M(3,5:5:5*(p+1)*(q+1)) = -theta*thickness/2*n1*NIJ;
                        F(dispIdxIJ) = F(dispIdxIJ) + M'* pload' * w3_vec(a);
                    end
                    
                end
            end
            
        end
        
        
    case 1 % lines force Fe = sum(M_GuassPoints'* P * det * w)  P = [fx,fy,fz] alone s = s0
        Weights = reshape(Ctrlpts(:,:,4), mAnalysis+1,nAnalysis+1);
        for kk = 1:size(ploadsMat,1)
            s0 = ploadsMat(kk,1);
            pload = ploadsMat(kk,3:5);
            
            for j = 1:ne_v
                [s2_vec, w_vec] = GaussInt(VEle(j),VEle(j+1),ng_v);
                uv = [ones(ng_v,1)*s0, s2_vec];
                
                Fe = zeros(5*(p+1)*(q+1),1);
                for k = 1:ng_v
                    [IdxU, IdxV, N, Ns, Nt, ~, ~, ~] = NrbSurDer2(mAnalysis,nAnalysis,p,q,uv(k,1),uv(k,2),UAnalysis,VAnalysis,Weights);
                    PIJ = Ctrlpts(IdxU+1,IdxV+1,1:3);
                    Ss = reshape(sum(PIJ.*Ns, [1,2]),1,3);
                    St = reshape(sum(PIJ.*Nt, [1,2]),1,3);
                    NIJ  = reshape(N,  1, []);
                    [J, I] = meshgrid(IdxV, IdxU);
                    Idx = sub2ind([mAnalysis+1, nAnalysis+1], I(:)+1, J(:)+1);
                    dispIdxIJ = reshape(((Idx-1)*5+(1:5))',1,[]);
                    
                    % local coordinate axises
                    v3 = cross(Ss,St);
                    v3 = v3 / norm(v3);
                    St_norm = norm(St);
                    Ss_norm = norm(Ss);
                    if(Ss_norm)
                        v1 = Ss/Ss_norm;
                        v2 = cross(v3, v1);
                    elseif(St_norm)
                        v2 = St/St_norm;
                        v1 = cross(v2, v3);
                    end
                    
                    l1 = v1(1); m1 = v1(2); n1 = v1(3);
                    l2 = v2(1); m2 = v2(2); n2 = v2(3);
                    
                    % Me
                    for a = 1:numel(w3_vec)
                        theta = s3_vec(a);
                        M = zeros(3,5*(p+1)*(q+1));
                        M(1,1:5:5*(p+1)*(q+1)) =  NIJ;
                        M(1,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*l2*NIJ;
                        M(1,5:5:5*(p+1)*(q+1)) =  theta*thickness/2*l1*NIJ;
                        M(2,2:5:5*(p+1)*(q+1)) =  NIJ;
                        M(2,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*m2*NIJ;
                        M(2,5:5:5*(p+1)*(q+1)) =  theta*thickness/2*m1*NIJ;
                        M(3,3:5:5*(p+1)*(q+1)) =  NIJ;
                        M(3,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*n2*NIJ;
                        M(3,5:5:5*(p+1)*(q+1)) = -theta*thickness/2*n1*NIJ;
                        
                        dl_dt = norm(St);  
                        Fe = Fe + M'* pload' * dl_dt * w_vec(k) * w3_vec(a);
                    end
                end
                
                F(dispIdxIJ) = F(dispIdxIJ) + Fe;
            end
            
        end
         
        
    case 2 %  lines force Fe = sum(M_GuassPoints'* P * det * w)  P = [fx,fy,fz] alone t = t0
        Weights = reshape(Ctrlpts(:,:,4), mAnalysis+1,nAnalysis+1);
        for kk = 1:size(ploadsMat,1)
            t0 = ploadsMat(kk,2);
            pload = ploadsMat(kk,3:5);
            
            for i = 1:ne_u
                [s1_vec, w_vec] = GaussInt(UEle(i),UEle(i+1),ng_u);
                uv = [s1_vec,ones(ng_u,1)*t0];
                
                Fe = zeros(5*(p+1)*(q+1),1);
                for k = 1:ng_u
                    [IdxU, IdxV, N, Ns, Nt, ~, ~, ~] = NrbSurDer2(mAnalysis,nAnalysis,p,q,uv(k,1),uv(k,2),UAnalysis,VAnalysis,Weights);
                    PIJ = Ctrlpts(IdxU+1,IdxV+1,1:3);
                    Ss = reshape(sum(PIJ.*Ns, [1,2]),1,3);
                    St = reshape(sum(PIJ.*Nt, [1,2]),1,3);
                    NIJ  = reshape(N,  1, []);
                    [J, I] = meshgrid(IdxV, IdxU);
                    Idx = sub2ind([mAnalysis+1, nAnalysis+1], I(:)+1, J(:)+1);
                    dispIdxIJ = reshape(((Idx-1)*5+(1:5))',1,[]);
                    
                    % local coordinate axises
                    v3 = cross(Ss,St);
                    v3 = v3 / norm(v3);
                    St_norm = norm(St);
                    Ss_norm = norm(Ss);
                    if(Ss_norm)
                        v1 = Ss/Ss_norm;
                        v2 = cross(v3, v1);
                    elseif(St_norm)
                        v2 = St/St_norm;
                        v1 = cross(v2, v3);
                    end
                    l1 = v1(1); m1 = v1(2); n1 = v1(3);
                    l2 = v2(1); m2 = v2(2); n2 = v2(3);
                    
                    % Me
                    for a = 1:numel(w3_vec)
                        theta = s3_vec(a);
                        M = zeros(3,5*(p+1)*(q+1));
                        M(1,1:5:5*(p+1)*(q+1)) =  NIJ;
                        M(1,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*l2*NIJ;
                        M(1,5:5:5*(p+1)*(q+1)) =  theta*thickness/2*l1*NIJ;
                        M(2,2:5:5*(p+1)*(q+1)) =  NIJ;
                        M(2,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*m2*NIJ;
                        M(2,5:5:5*(p+1)*(q+1)) =  theta*thickness/2*m1*NIJ;
                        M(3,3:5:5*(p+1)*(q+1)) =  NIJ;
                        M(3,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*n2*NIJ;
                        M(3,5:5:5*(p+1)*(q+1)) = -theta*thickness/2*n1*NIJ;
                        
                        dl_ds = norm(Ss);
                        Fe = Fe + M'* pload' * dl_ds * w_vec(k) * w3_vec(a);
                    end
                end
                
                F(dispIdxIJ) = F(dispIdxIJ) + Fe;
            end
        end
        
        
    case 3 %  lines force Fe = sum(M_GuassPoints'* P * det * w)  P = [fx,fy,fz] alone s = s0 and t = t0
        Weights = reshape(Ctrlpts(:,:,4), mAnalysis+1,nAnalysis+1);
        for kk = 1:size(ploadsMat,1)
            s0 = ploadsMat(kk,1);
            t0 = ploadsMat(kk,2);
            pload = ploadsMat(kk,3:5);
            
            for j = 1:ne_v
                [s2_vec, w_vec] = GaussInt(VEle(j),VEle(j+1),ng_v);
                uv = [ones(ng_v,1)*s0, s2_vec];
                
                Fe = zeros(5*(p+1)*(q+1),1);
                for k = 1:ng_v
                    [IdxU, IdxV, N, Ns, Nt, ~, ~, ~] = NrbSurDer2(mAnalysis,nAnalysis,p,q,uv(k,1),uv(k,2),UAnalysis,VAnalysis,Weights);
                    PIJ = Ctrlpts(IdxU+1,IdxV+1,1:3);
                    Ss = reshape(sum(PIJ.*Ns, [1,2]),1,3);
                    St = reshape(sum(PIJ.*Nt, [1,2]),1,3);
                    NIJ  = reshape(N,  1, []);
                    [J, I] = meshgrid(IdxV, IdxU);
                    Idx = sub2ind([mAnalysis+1, nAnalysis+1], I(:)+1, J(:)+1);
                    dispIdxIJ = reshape(((Idx-1)*5+(1:5))',1,[]);
                    
                    % local coordinate axises
                    v3 = cross(Ss,St);
                    v3 = v3 / norm(v3);
                    St_norm = norm(St);
                    Ss_norm = norm(Ss);
                    if(Ss_norm)
                        v1 = Ss/Ss_norm;
                        v2 = cross(v3, v1);
                    elseif(St_norm)
                        v2 = St/St_norm;
                        v1 = cross(v2, v3);
                    end
                    l1 = v1(1); m1 = v1(2); n1 = v1(3);
                    l2 = v2(1); m2 = v2(2); n2 = v2(3);
                    
                    % Me
                    for a = 1:numel(w3_vec)
                        theta = s3_vec(a);
                        M = zeros(3,5*(p+1)*(q+1));
                        M(1,1:5:5*(p+1)*(q+1)) =  NIJ;
                        M(1,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*l2*NIJ;
                        M(1,5:5:5*(p+1)*(q+1)) =  theta*thickness/2*l1*NIJ;
                        M(2,2:5:5*(p+1)*(q+1)) =  NIJ;
                        M(2,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*m2*NIJ;
                        M(2,5:5:5*(p+1)*(q+1)) =  theta*thickness/2*m1*NIJ;
                        M(3,3:5:5*(p+1)*(q+1)) =  NIJ;
                        M(3,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*n2*NIJ;
                        M(3,5:5:5*(p+1)*(q+1)) = -theta*thickness/2*n1*NIJ;
                        
                        dl_dt = norm(St);
                        Fe = Fe + M'* pload' * dl_dt * w_vec(k) * w3_vec(a);
                    end
                end
                
                F(dispIdxIJ) = F(dispIdxIJ) + Fe;
            end
            
            for i = 1:ne_u
                [s1_vec, w_vec] = GaussInt(UEle(i),UEle(i+1),ng_u);
                uv = [s1_vec,ones(ng_u,1)*t0];
                
                Fe = zeros(5*(p+1)*(q+1),1);
                for k = 1:ng_u
                    [IdxU, IdxV, N, Ns, Nt, ~, ~, ~] = NrbSurDer2(mAnalysis,nAnalysis,p,q,uv(k,1),uv(k,2),UAnalysis,VAnalysis,Weights);
                    PIJ = Ctrlpts(IdxU+1,IdxV+1,1:3);
                    Ss = reshape(sum(PIJ.*Ns, [1,2]),1,3);
                    St = reshape(sum(PIJ.*Nt, [1,2]),1,3);
                    NIJ  = reshape(N,  1, []);
                    [J, I] = meshgrid(IdxV, IdxU);
                    Idx = sub2ind([mAnalysis+1, nAnalysis+1], I(:)+1, J(:)+1);
                    dispIdxIJ = reshape(((Idx-1)*5+(1:5))',1,[]);
                    
                    % local coordinate axises
                    v3 = cross(Ss,St);
                    v3 = v3 / norm(v3);
                    St_norm = norm(St);
                    Ss_norm = norm(Ss);
                    if(Ss_norm)
                        v1 = Ss/Ss_norm;
                        v2 = cross(v3, v1);
                    elseif(St_norm)
                        v2 = St/St_norm;
                        v1 = cross(v2, v3);
                    end
                    l1 = v1(1); m1 = v1(2); n1 = v1(3);
                    l2 = v2(1); m2 = v2(2); n2 = v2(3);
                    
                    % Me
                    for a = 1:numel(w3_vec)
                        theta = s3_vec(a);
                        M = zeros(3,5*(p+1)*(q+1));
                        M(1,1:5:5*(p+1)*(q+1)) =  NIJ;
                        M(1,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*l2*NIJ;
                        M(1,5:5:5*(p+1)*(q+1)) =  theta*thickness/2*l1*NIJ;
                        M(2,2:5:5*(p+1)*(q+1)) =  NIJ;
                        M(2,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*m2*NIJ;
                        M(2,5:5:5*(p+1)*(q+1)) =  theta*thickness/2*m1*NIJ;
                        M(3,3:5:5*(p+1)*(q+1)) =  NIJ;
                        M(3,4:5:5*(p+1)*(q+1)) = -theta*thickness/2*n2*NIJ;
                        M(3,5:5:5*(p+1)*(q+1)) = -theta*thickness/2*n1*NIJ;
                        
                        dl_ds = norm(Ss);
                        Fe = Fe + M'* pload' * dl_ds * w_vec(k) * w3_vec(a);
                    end
                end
                
                F(dispIdxIJ) = F(dispIdxIJ) + Fe;
            end
            
        end
        
end

end