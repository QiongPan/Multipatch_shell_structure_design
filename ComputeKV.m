/**
 * @file ComputeKV.m
 * @author Qiong Pan (PQ2019@mail.ustc.edu.cn)
 * @brief  compute the stiffness and volume on guass points (elementwise) for shell analysis.
 * Orginal work written by Qiong Pan in Matlab.
 * @date 2024-01-23
 */

function [Ke_GuassPt_ECell, Ve_GuassPt_All,Weight_GuassPt_All] = ComputeKV(mAnalysis,nAnalysis,p,q,UAnalysis,VAnalysis,UEle,VEle,Ctrlpts, ne_u, ne_v)

%%%%%% stiffness matrix for each patch, computed elementwise on Guass integration points %%%%%%
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
%   ne_u       - number of element(knot spans) in u-direction    (U-direction)
%   ne_v       - number of element(knot spans) in v-direction    (V-direction)
% output:
%   Ke_GuassPt_ECell    - stiffness matrix in Gauss integration points of elements
%                       - a cell with dim = (ne_u*ne_v) * 1
%                       - each cell element stores a matrix with dim = (ng_u*ng_v) * (5*(p+1)*(q+1))^2)
%                       - each row of the matrix is an expansion of the stiffness matrix at a Gaussian integration point
%   Ve_GuassPt_All      - volume for Gauss integration
%                       - a matrix with dim = (ne_u*ne_v) * (ng_u*ng_v)
%   Weight_GuassPt_All  - weights for Gauss integration
%                       - a matrix with dim = (ne_u*ne_v) * (ng_u*ng_v)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K = \int\limits_{\Omega} {B'*D*B} d{\Omega}
%   = \int\limits_{0}^1 \int\limits_{0}^1 \int\limits_{-1}^1 {{B_k}'*D*{B_k}} * |det(J)|  d{\zeta} dt ds
%   = \sum\limits_e {\int\limits_{\Omega_e} {B'*D*B} d{\Omega_e}}                   -- elementwise
%   = \sum\limits_e {\int\limits_{UEle(i)}^{UEle(i)} \int\limits_{VEle(j)}^{VEle(j)} \int\limits_{-1}^1 B'*D*B * |det(J)| d{\zeta} dt ds
%   = \sum\limits_e {\sum\limits_k^{N_G} {{B_k}'*D*{B_k} * |det(J_k)| W_k }  }      -- Guass integration
%
% B = T * H * \Gamma * R
%   - T: transformation matrix between local and global strain.
%        dim = 6 * 6
%   - H: a constant matrix: ({\epsilon}_{xx}, {\epsilon}_{yy}, ... {\gamma}_{xz}) = H * (du/dx, du/dy, du/dz, dv/dx, ... dw/dz)
%        dim = 6 * 9
%   - \Gamma: a matrix that transforms the partial derivatives to those w.r.t parameters (s, t, \zeta)
%        dim = 9 * 9
%   - R: a matrix that correlates these partial derivatives w.r.t.(s, t, \zeta) with displacement variables.
%        dim = 9 * (5*(p+1)*(q+1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  To distinguish the surface symbol 'S', the knot vectors are represented by U and V.
%  Thus, u and s represent the same parameter direction, the relationship between t and v is similar.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



global poisson thickness ng_u ng_v ng_th;
% poisson   - Poisson's ratio
% thickness - the thickness 'h' of the shell
% ng_u      - number of Gauss integration point along u-direction
% ng_v      - number of Gauss integration point along v-direction
% ng_th     - number of Gauss integration point along thickness direction

Weights             = reshape(Ctrlpts(:,:,4),mAnalysis+1,nAnalysis+1);  % weights of NURBS for the mid-surface
Ke_GuassPt_ECell    = cell(ne_u*ne_v, 1);                               % ele_idx = (j-1)*ne_u + i; gp = [s1,t1; s2 t1; ...; s1 t2,...]; {Ke(ng_u*ng_v,(5*(p+1)*(q+1))^2)}
Ve_GuassPt_All      = zeros(ne_u*ne_v, ng_u*ng_v);
Weight_GuassPt_All  = zeros(ne_u*ne_v, ng_u*ng_v);

[H, D]  = ConstMat(1, poisson);             % Constant matrix
[s3_vec,w3_vec] = GaussInt(-1,1,ng_th);     % Guass intrgration points along thickness

for j = 1:ne_v
    [s2_vec, w2_vec] = GaussInt(VEle(j),VEle(j+1),ng_v);        % Guass intrgration points along v-direction
    for i = 1:ne_u
        [s1_vec, w1_vec] = GaussInt(UEle(i),UEle(i+1),ng_u);    % Guass intrgration points along u-direction
        [tmp1, tmp2] = meshgrid(s1_vec,s2_vec);
        uv = [reshape(tmp1,[],1),reshape(tmp2,[],1)];           % ordered vector of Guass intrgration points in parametric space [0,1]^2: gp = [s1,t1; s2 t1; ...; s1 t2,...];
        [tmp1, tmp2] = meshgrid(w1_vec,w2_vec);
        w = reshape(tmp1.*tmp2,1,[]);                           % weight vector of Guass intrgration points
        ele_idx = (j-1)*ne_u+i;
        Ke_GuassPt = zeros(ng_u*ng_v, (5*(p+1)*(q+1))^2);       % a matrix, each row contains the expansion of the stiffness matrix at a Gaussian integration point
        
        for k = 1:numel(w)
            % KeGP & VeGP
            VeGP = 0;
            KeGP = zeros(5*(p+1)*(q+1));
            
            [IdxU, IdxV, N, Ns, Nt, Nst, Nss, Ntt] = NrbSurDer2(mAnalysis,nAnalysis,p,q,uv(k,1),uv(k,2),UAnalysis,VAnalysis,Weights);
            PIJ = Ctrlpts(IdxU+1,IdxV+1,1:3);
            Ss = reshape(sum(PIJ.*Ns, [1,2]),1,3);
            St = reshape(sum(PIJ.*Nt, [1,2]),1,3);
            Sss = reshape(sum(PIJ.*Nss, [1,2]),1,3);
            Stt = reshape(sum(PIJ.*Ntt, [1,2]),1,3);
            Sst = reshape(sum(PIJ.*Nst, [1,2]),1,3);
            
            %% local coordinate axises
            v3 = cross(Ss,St);
            v3_norm = norm(v3);
            v3_s = ((cross(Sss,St)+cross(Ss,Sst))*(v3_norm^2) - v3*dot(v3,cross(Sss,St)+cross(Ss,Sst)))/(v3_norm^3);
            v3_t = ((cross(Ss,Stt)+cross(Sst,St))*(v3_norm^2) - v3*dot(v3,cross(Sst,St)+cross(Ss,Stt)))/(v3_norm^3);
            v3 = v3 / v3_norm;
            St_norm = norm(St);
            Ss_norm = norm(Ss);
            if(Ss_norm)
                v1 = Ss/Ss_norm;
                v1_s = (Sss - v1*dot(Ss,Sss)/Ss_norm) / Ss_norm;
                v1_t = (Sst - v1*dot(Ss,Sst)/Ss_norm) / Ss_norm;
                v2 = cross(v3, v1);
                v2_s = cross(v3_s, v1) + cross(v3, v1_s);
                v2_t = cross(v3_t, v1) + cross(v3, v1_t);
            elseif(St_norm)
                v2 = St/St_norm;
                v2_s = (Sst - v2*dot(St,Sst)/St_norm) / St_norm;
                v2_t = (Stt - v2*dot(St,Stt)/St_norm) / St_norm;
                v1 = cross(v2, v3);
                v1_s = cross(v2_s, v3) + cross(v2, v3_s);
                v1_t = cross(v2_t, v3) + cross(v2, v3_t);
            end
            
            
            %%
            l1 = v1(1); m1 = v1(2); n1 = v1(3);
            l2 = v2(1); m2 = v2(2); n2 = v2(3);
            l3 = v3(1); m3 = v3(2); n3 = v3(3);
            l1_s = v1_s(1); m1_s = v1_s(2); n1_s = v1_s(3);
            l1_t = v1_t(1); m1_t = v1_t(2); n1_t = v1_t(3);
            l2_s = v2_s(1); m2_s = v2_s(2); n2_s = v2_s(3);
            l2_t = v2_t(1); m2_t = v2_t(2); n2_t = v2_t(3);
            
            %% matrix B = THGR
            T = [ l1^2       m1^2       n1^2       l1*m1       m1*n1       n1*l1
                l2^2       m2^2       n2^2       l2*m2       m2*n2       n2*l2
                l3^2       m3^2       n3^2       l3*m3       m3*n3       n3*l3
                2*l1*l2    2*m1*m2    2*n1*n2    l1*m2+l2*m1 m1*n2+m2*n1 n1*l2+n2*l1
                2*l2*l3    2*m2*m3    2*n2*n3    l2*m3+l3*m2 m2*n3+m3*n2 n2*l3+n3*l2
                2*l3*l1    2*m3*m1    2*n3*n1    l3*m1+l1*m3 m3*n1+m1*n3 n3*l1+n1*l3];
            
            NIJ  = reshape(N,  1, []);
            NsIJ = reshape(Ns, 1, []);
            NtIJ = reshape(Nt, 1, []);
            
            %% theta-direction integration;
            for a = 1:numel(w3_vec)
                theta = s3_vec(a);
                JS = [Ss; St; zeros(1,3)] + thickness/2*[theta*v3_s; theta*v3_t; v3];
                invJS = inv(JS);
                detJS = det(JS);
                
                Gamma = blkdiag(invJS,invJS,invJS);
                
                %% compute R
                R = zeros(9,5*(p+1)*(q+1));
                R(1, 1:5:end) =  NsIJ;
                R(1, 4:5:end) = -theta*thickness/2*(l2_s*NIJ+l2*NsIJ);
                R(1, 5:5:end) =  theta*thickness/2*(l1_s*NIJ+l1*NsIJ);
                R(2, 1:5:end) =  NtIJ;
                R(2, 4:5:end) = -theta*thickness/2*(l2_t*NIJ+l2*NtIJ);
                R(2, 5:5:end) =  theta*thickness/2*(l1_t*NIJ+l1*NtIJ);
                R(3, 4:5:end) = -thickness/2*(l2*NIJ);
                R(3, 5:5:end) =  thickness/2*(l1*NIJ);
                R(4, 2:5:end) =  NsIJ;
                R(4, 4:5:end) = -theta*thickness/2*(m2_s*NIJ+m2*NsIJ);
                R(4, 5:5:end) =  theta*thickness/2*(m1_s*NIJ+m1*NsIJ);
                R(5, 2:5:end) =  NtIJ;
                R(5, 4:5:end) = -theta*thickness/2*(m2_t*NIJ+m2*NtIJ);
                R(5, 5:5:end) =  theta*thickness/2*(m1_t*NIJ+m1*NtIJ);
                R(6, 4:5:end) = -thickness/2*(m2*NIJ);
                R(6, 5:5:end) =  thickness/2*(m1*NIJ);
                R(7, 3:5:end) =  NsIJ;
                R(7, 4:5:end) = -theta*thickness/2*(n2_s*NIJ+n2*NsIJ);
                R(7, 5:5:end) =  theta*thickness/2*(n1_s*NIJ+n1*NsIJ);
                R(8, 3:5:end) =  NtIJ;
                R(8, 4:5:end) = -theta*thickness/2*(n2_t*NIJ+n2*NtIJ);
                R(8, 5:5:end) =  theta*thickness/2*(n1_t*NIJ+n1*NtIJ);
                R(9, 4:5:end) = -thickness/2*(n2*NIJ);
                R(9, 5:5:end) =  thickness/2*(n1*NIJ);
                
                %%
                B = T * H * Gamma * R;
                KeGP = KeGP + B' * D * B * detJS * w3_vec(a);
                VeGP = VeGP + detJS * w3_vec(a);
            end
            
            Ke_GuassPt(k,:) = reshape(KeGP,[],1);
            Ve_GuassPt_All(ele_idx,k) = VeGP;
        end
        Ke_GuassPt_ECell(ele_idx) = {Ke_GuassPt};
        Weight_GuassPt_All(ele_idx,:) = w;
        
    end
end

end
