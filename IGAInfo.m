/**
 * @file IGAInfo.m
 * @author Qiong Pan (PQ2019@mail.ustc.edu.cn)
 * @brief analysis and sensitivity analysis for a multi-patch shell structure.
 * Orginal work written by Qiong Pan in Matlab.
 * @date 2024-01-23
 */

function [compl,DersCompl_Vec,vol,DersVol_Vec,GlobalVol] = IGAInfo(DOF_Analysis,DOF_Design,num_pat,DegreeMat,num_ele_Mat,DensCtrlpts,NonZeroBasis_GuassPt_PCell,DNonZeroIdx_PCell, ...
    ANonZeroIdx_PCell,StiffnessVolMat_PCell,F_Global,fixed_DOF,myFilter,EleInfoMat,EleIdx_PCell,Neighbors_ECell)

%%%%%% compute compliance, volume and their derivatives w.r.t. density coefficiences %%%%%%
% input:
%   DOF_Analysis                - degrees of freedom in analysis model
%                               - number of control points in analysis model * 5
%   DOF_Design                  - degrees of freedom in design model
%                               - number of control points/coefficients in design model
%   num_pat                     - number of matching patches
%   DegreeMat                   - DegreeMat(k,:) = (p_k, q_k), degree matrix for patches
%                               - dim = num_pat * 2;
%   num_ele_Mat                 - num_ele_Mat(k,:) = (ne_u_k, ne_v_k)
%                               - dim = num_pat * 2;
%   DensCtrlpts                 - density coefficients, i.e. design variables
%                               - dim = DOF_Design * 1;
%   NonZeroBasis_GuassPt_PCell 	- nonzero basis functions in design model at Gauss points of elements of patches
%                               - a cell with dim = num_pat * 1 (patch related)
%                               - each patch related cell element stores a cell with dim = (ne_u*ne_v) * 1 (element related)
%                               - each element related cell stores a =(ng_u*ng_t) * ((p+1)*(q+1)) matrix
%                               - each row of the matrix records the nonzero basis functions in design model at a Gauss point
%   DNonZeroIdx_PCell           - indices of nonzero basis functions in design model at Gauss points of elements of patches
%                               - a cell with dim = num_pat * 1 (patch related)
%                               - each patch related cell element stores a cell with dim = (ne_u*ne_v) * 1 (element related)
%                               - each element related cell stores a =(ng_u*ng_t) * ((p+1)*(q+1)) matrix
%                               - each row of the matrix records the indicies of nonzero basis functions at a Gauss point
%   ANonZeroIdx_PCell           - indices of nonzero basis functions in analysis model in design model at Gauss points of elements of patches
%                               - a cell with dim = num_pat * 1
%                               - each cell element stores a matrix with dim = (ne_u*ne_v) * ((p+1)*(q+1))
%                               - each row of the matrix records the indicies of nonzero basis functions in analysis model at a Gauss point
%                               - dim = DOF_Design * 1;
%   StiffnessVolMat_PCell       - a patch related cell that stores the stiffness matrix, volume and weights at guass points
%                               - StiffnessVolWeight_PCell(k,:) = {Ke_GuassPt_ECell, Ve_GuassPt_All, Weight_GuassPt_All}
%                               - Ke_GuassPt_ECell  : stiffness matrix in Gauss integration points of elements (element ralated cell)
%                               - Ve_GuassPt_All    : volume for Gauss integration (element ralated matrix)
%                               - Weight_GuassPt_All: weights for Gauss integration (element ralated matrix)
%   F_Global                    - Forcr vector F (K*U = F)
%                               - dim = DOF_Analysis * 1;
%   fixed_DOF                   - fixed dofs in analysis model
%   myFilter                    - 'NULL': total volume constraint with no filter
%                               - 'LocalVolumeConstraint': |Pi-Pj| <= Radius = filter_adius * elements_edge_edgelength_average
%   EleInfoMat                  - element information matrix with EleInfoMat(ele,:) = (ele_pat, ele_prow, ele_pcol)
%                               - each row illustrates the patch to which the element belongs, the row index and column index of the element in that patch
%                               - dim = (number of elements in total) * 3
%                               - for Local Volume Constraint only
%   EleIdx_PCell                - the global index for elements in each patch. (patch related cell)
%                               - dim = (number of elements in total) * 1
%                               - each cell element stores a ne_u_k * ne_v_k matrix
%                               - for Local Volume Constraint only
%   Neighbors_ECell             - the indicies of neighbor elements for each element in patches(element related matrix)
%                               - dim = (number of elements in total) * 1
%                               - each cell element stores a cell
%                               - for Local Volume Constraint only
% output:
%   compl                       - compliance of the optimized shell structure
%   DersCompl_Vec               - the partial derivatives of compl w.r.t. density coefficients
%   vol                         - total volume of the optimized shell structure for volume fraction constraint
%                               - local volume fraction(the maximum value in a p-norm sense)
%   DersVol_Vec                 - the partial derivatives of vol w.r.t. density coefficients
%   GlobalVol                   - total volume of the optimized shell structure
%                               - used in local volume constraint to record the total material volume

Dense_GuassPt_All_PCell = ComputeDen(num_pat,num_ele_Mat,DensCtrlpts,NonZeroBasis_GuassPt_PCell,DNonZeroIdx_PCell);
K_Global                = ComputeK(DegreeMat,DOF_Analysis,num_pat,num_ele_Mat,StiffnessVolMat_PCell,Dense_GuassPt_All_PCell,ANonZeroIdx_PCell);
UCtrlpts                = SolveU(DOF_Analysis, K_Global, F_Global,fixed_DOF);
clear DOF_Analysis F_Global DispC1_Mat F_Global fixed_DOF;

GlobalVol = 0;
switch myFilter
    case 'NULL'
        [compl,DersCompl_Vec,vol,DersVol_Vec] = ComputeDerCV(num_pat,DOF_Design,DegreeMat,num_ele_Mat,UCtrlpts,Dense_GuassPt_All_PCell,StiffnessVolMat_PCell,...
            NonZeroBasis_GuassPt_PCell,DNonZeroIdx_PCell,ANonZeroIdx_PCell);
    case 'LocalVolumeConstraint'
        [compl,DersCompl_Vec,vol,DersVol_Vec,GlobalVol] = ComputeDerCV_LVC(num_pat,DOF_Design,DegreeMat,num_ele_Mat,UCtrlpts,Dense_GuassPt_All_PCell,StiffnessVolMat_PCell,...
            NonZeroBasis_GuassPt_PCell,DNonZeroIdx_PCell,ANonZeroIdx_PCell,EleInfoMat, EleIdx_PCell, Neighbors_ECell);
end

end

%% density for guass points in elenments
function Dense_GuassPt_All_PCell  = ComputeDen(num_pat,num_ele_Mat,DensCtrlpts,NonZeroBasis_GuassPt_PCell,DNonZeroIdx_PCell)
Dense_GuassPt_All_PCell = cell(num_pat,1);
global ng_u ng_v;
for k = 1:num_pat
    ne_u = num_ele_Mat(k,1);
    ne_v = num_ele_Mat(k,2);
    NonZeroBasis_GuassPt_ECell = NonZeroBasis_GuassPt_PCell{k};
    DNonZeroIdx_ECell = DNonZeroIdx_PCell{k};
    Dense_GuassPt_All = zeros(ne_u*ne_v,ng_u*ng_v);
    for j = 1:ne_v
        for i = 1:ne_u
            ij = (j-1)*ne_u+i;
            for kgp = 1:ng_u*ng_v
                dense = dot(NonZeroBasis_GuassPt_ECell{ij}(kgp,:), DensCtrlpts(DNonZeroIdx_ECell{ij}(kgp,:)));
                Dense_GuassPt_All(ij,kgp) = dense;
            end
        end
    end
    Dense_GuassPt_All_PCell(k) = {Dense_GuassPt_All};
end

end

%% stiffness
function K_Global = ComputeK(DegreeMat,DOF_Analysis,num_pat,num_ele_Mat,StiffnessVolMat_PCell,Dense_All_PCell,ANonZeroIdx_PCell)
global penal E0 Emin beta ng_u ng_v;
ngst = ng_u*ng_v;
iH= [];
jH= [];
sH= [];
for k = 1:num_pat
    ne_u = num_ele_Mat(k,1);
    ne_v = num_ele_Mat(k,2);
    Ke_ECell = StiffnessVolMat_PCell{k,1};
    Weight_All = StiffnessVolMat_PCell{k,3};
    Dense_All = Dense_All_PCell{k};
    ANonZeroIdx = ANonZeroIdx_PCell{k};
    
    itr = 0;
    p = DegreeMat(k,1);
    q = DegreeMat(k,2);
    pq = (5*(p+1)*(q+1));
    num_ =  ne_u*ne_v;
    RowMat = zeros(pq^2, num_);
    ColMat = zeros(pq^2, num_);
    ValMat = zeros(pq^2, num_);
    
    for j = 1:ne_v
        for i = 1:ne_u
            ij = (j-1)*ne_u+i;
            IdxIJ = ANonZeroIdx(ij,:);
            KIdx  = reshape((IdxIJ-1)*5+(1:5)',[],1);
            [Y,X] = meshgrid(KIdx,KIdx);
            
            itr = itr + 1;
            RowMat(:,itr) = reshape(X,pq^2,1);
            ColMat(:,itr) = reshape(Y,pq^2,1);
            Ke = Ke_ECell{ij};
            
            for kgp = 1:ngst
                Dense = Dense_All(ij,kgp);
                HDense = Heaviside(Dense,beta,0.5);
                E = Emin + HDense^penal*(E0-Emin);
                ValMat(:,itr) = ValMat(:,itr) + Ke(kgp,:)'*E*Weight_All(ij,kgp);
            end
        end
    end
    
    iH = [iH;reshape(RowMat,[],1)];
    jH = [jH;reshape(ColMat,[],1)];
    sH = [sH;reshape(ValMat,[],1)];
end

K_Global = sparse(iH,jH,sH,5*DOF_Analysis,5*DOF_Analysis);

end

%% Control Points of Displacement
function UCtrlpts = SolveU(DOF_Analysis, K_Global, F_Global,fixed_DOF)
UCtrlpts = zeros(5*DOF_Analysis,1);	% control coefficients of diaplacement U: u_ij
K_Global = (K_Global+K_Global')/2;

fixed = reshape(((fixed_DOF-1)*5+(1:5)'),[],1);
freedofs = setdiff((1:5*DOF_Analysis)', fixed);
K_      = K_Global(freedofs,freedofs);
F_      = F_Global(freedofs); % COND = condest(K_)
Ufree = K_ \ F_;
UCtrlpts(freedofs) = Ufree;

% Ufree = minres(K_, F_);
% Ufree = minres(K_, F_, 1e-6, 500);

% % QR division
% [CC,RR,PP] = qr(K_,F_);
% Ufree = PP*(RR\CC);

% Ufree = lsqlin(K_,F_,[],[]);
% UCtrlpts(setdiff([1:5*DOF_Analysis]', fixed)) = Ufree;
end

%% compliance,volume and their derivatives with Density-coefficients filter
function [compl,DersCompl_Vec,vol,DersVol_Vec] = ComputeDerCV(num_pat,DOF_Design,DegreeMat,num_ele_Mat,UCtrlpts,Dense_All_PCell,StiffnessVolMat_PCell,NonZeroBasis_GuassPt_PCell,DNonZeroIdx_PCell,ANonZeroIdx_PCell)
global beta penal E0 Emin ng_u ng_v;
compl = 0.0;
vol   = 0.0;
DersCompl_Vec = zeros(DOF_Design,1);
DersVol_Vec = zeros(DOF_Design,1);

for k = 1:num_pat
    p = DegreeMat(k,1);
    q = DegreeMat(k,2);
    Dense_All = Dense_All_PCell{k};
    Ke_ECell = StiffnessVolMat_PCell{k,1};
    Ve_All = StiffnessVolMat_PCell{k,2};
    Weight_All = StiffnessVolMat_PCell{k,3};
    NonZeroBasis_GuassPt_ECell = NonZeroBasis_GuassPt_PCell{k};
    DNonZeroIdx_ECell = DNonZeroIdx_PCell{k};
    ANonZeroIdx = ANonZeroIdx_PCell{k};
    
    ne_u = num_ele_Mat(k,1);
    ne_v = num_ele_Mat(k,2);
    for j = 1:ne_v
        for i = 1:ne_u
            ij = (j-1)*ne_u + i;
            AIdxIJ = ANonZeroIdx(ij,:);
            UIdx  = reshape((AIdxIJ-1)*5+(1:5)',[],1);
            Ue = UCtrlpts(UIdx);
            
            for kgp = 1:ng_u*ng_v
                Ve = Ve_All(ij,kgp);
                Dense = Dense_All(ij,kgp);
                HDense = Heaviside(Dense,beta,0.5);
                DerHDense = DerHeaviside(Dense,beta,0.5);
                
                Ke0   = reshape(Ke_ECell{ij}(kgp,:),5*(p+1)*(q+1),5*(p+1)*(q+1))*Weight_All(ij,kgp);
                E = Emin + HDense^penal*(E0-Emin);
                UKU = (Ue' * Ke0 * Ue);
                compl = compl + E * UKU;
                vol = vol + Ve * HDense * Weight_All(ij,kgp);
                
                
                DNodeIdx = DNonZeroIdx_ECell{ij}(kgp,:)';
                NIJ = NonZeroBasis_GuassPt_ECell{ij}(kgp,:)';
                DersCompl_Vec(DNodeIdx) = DersCompl_Vec(DNodeIdx) - UKU*(E0-Emin) * (penal * HDense^(penal-1) * DerHDense * NIJ);
                DersVol_Vec(DNodeIdx) = DersVol_Vec(DNodeIdx) + Ve * Weight_All(ij,kgp) * DerHDense * NIJ;
            end
            
        end
    end
end

end

%% compliance,volume and their derivatives with Local Volume Constraint
function [compl,DersCompl_Vec,LocalVol,DersLocalVol, GlobalVol] = ComputeDerCV_LVC(num_pat,DOF_Design,DegreeMat,num_ele_Mat,UCtrlpts,Dense_All_PCell,StiffnessVolMat_PCell,NonZeroBasis_GuassPt_PCell,DNonZeroIdx_PCell,ANonZeroIdx_PCell,EleInfoMat, EleIdx_PCell, Neighbors_ECell)
% EleInfoMat = [];                  % Record Element Information: (ele_pat, ele_prow, ele_pcol).
% Neighbors_ECell = [];             % Record the indicies of neighbor-elements for each element.
% EleIdx_PCell = cell(num_pat,1);	% Record the global ele-idx for elements of each patch.
global beta penal E0 Emin ng_u ng_v;
ngst = ng_u * ng_v;
compl       = 0;
LocalVol    = 0;
GlobalVol   = 0;
DersCompl_Vec = zeros(DOF_Design,1);
DersLocalVol = zeros(DOF_Design,1);

% Denseity of elements(guass points) and derivatives
num_ele = size(Neighbors_ECell,1);
HDense_All = zeros(num_ele,ngst);
DerHDense_All = HDense_All;
HDense_GP_Vec = zeros(1,ngst);
DerHDense_GP_Vec = zeros(1,ngst);
for k = 1:num_pat
    EleIdx_Mat = EleIdx_PCell{k};
    Dense_All = Dense_All_PCell{k};
    Ve_All = StiffnessVolMat_PCell{k,2};
    Weight_All = StiffnessVolMat_PCell{k,3};
    GlobalVol = GlobalVol + sum(Ve_All.*Dense_All.*Weight_All,'all');
    ne_u = num_ele_Mat(k,1);
    ne_v = num_ele_Mat(k,2);
    for j = 1:ne_v
        for i = 1:ne_u
            ij = (j-1)*ne_u + i;
            for kgp = 1:ngst
                Dense = Dense_All(ij,kgp);
                HDense_GP_Vec(kgp) = Heaviside(Dense,beta,0.5);
                DerHDense_GP_Vec(kgp) = DerHeaviside(Dense,beta,0.5);
            end
            
            HDense_All(EleIdx_Mat(i,j),:) = HDense_GP_Vec;
            DerHDense_All(EleIdx_Mat(i,j),:) = DerHDense_GP_Vec;
        end
    end
end

% local volume of elements
LocalVe_All = zeros(num_ele,1);     % g_e = g(rho_e) = Simgma_e ((H(rho_e))*Ve) / Simgma_e Ve
SumLocalVe0_All = zeros(num_ele,1); % Simgma_e(Ve)
for ele_idx = 1:num_ele
    LocalVe = 0;
    SumLocalVe0 = 0;
    Neighbors_idx = Neighbors_ECell{ele_idx};
    for nidx = 1:numel(Neighbors_idx)  % global idx for this element
        pat_idx = EleInfoMat(Neighbors_idx(nidx),1);
        row_idx = EleInfoMat(Neighbors_idx(nidx),2);
        col_idx = EleInfoMat(Neighbors_idx(nidx),3);
        ne_u = num_ele_Mat(pat_idx,1);
        ij = (col_idx-1)*ne_u + row_idx; % local idx within the patch: pat_idx, for this element
        Ve_All = StiffnessVolMat_PCell{pat_idx,2};
        Weight_All = StiffnessVolMat_PCell{pat_idx,3};
        LocalVe = LocalVe + sum(HDense_All(Neighbors_idx(nidx),:).*Ve_All(ij,:).*Weight_All(ij,:),'all');
        SumLocalVe0 = SumLocalVe0 + sum(Ve_All(ij,:).*Weight_All(ij,:),'all');
    end
    LocalVe_All(ele_idx) = LocalVe / SumLocalVe0;
    SumLocalVe0_All(ele_idx) = SumLocalVe0;
end

% Compliance and derivatives w.r.t. \rho_{ij}
for k = 1:num_pat
    p = DegreeMat(k,1);
    q = DegreeMat(k,2);
    EleIdx_Mat = EleIdx_PCell{k};
    Ke_ECell = StiffnessVolMat_PCell{k,1};
    Weight_All = StiffnessVolMat_PCell{k,3};
    NonZeroBasis_GuassPt_ECell = NonZeroBasis_GuassPt_PCell{k};
    DNonZeroIdx_ECell = DNonZeroIdx_PCell{k};
    ANonZeroIdx = ANonZeroIdx_PCell{k};
    
    ne_u = num_ele_Mat(k,1);
    ne_v = num_ele_Mat(k,2);
    for j = 1:ne_v
        for i = 1:ne_u
            ij = (j-1)*ne_u + i;
            ele_idx = EleIdx_Mat(i,j); % global idx
            AIdxIJ = ANonZeroIdx(ij,:);
            UIdx  = reshape((AIdxIJ-1)*5+(1:5)',[],1);
            Ue = UCtrlpts(UIdx);
            
            for kgp = 1:ngst
                HDense = HDense_All(ele_idx,kgp);
                DerHDense = DerHDense_All(ele_idx,kgp);
                
                Ke0   = reshape(Ke_ECell{ij}(kgp,:),5*(p+1)*(q+1),5*(p+1)*(q+1))*Weight_All(ij,kgp);
                E = Emin + HDense^penal*(E0-Emin);
                UKU = (Ue' * Ke0 * Ue);
                compl = compl + E * UKU;
                
                % derivative of Compliance w.r.t. \rho_{ij}
                DNodeIdx = DNonZeroIdx_ECell{ij}(kgp,:)';
                NIJ = NonZeroBasis_GuassPt_ECell{ij}(kgp,:)';
                DersCompl_Vec(DNodeIdx) = DersCompl_Vec(DNodeIdx) - UKU*(E0-Emin) * (penal * HDense^(penal-1) * DerHDense * NIJ);
            end
        end
    end
end

% local volume and derivatives  w.r.t. \rho_{ij}
pow = 64;  % power
SumLVe_pow = sum(LocalVe_All.^pow,'all');   % Simgma_e (g_e)^p
LocalVol = (SumLVe_pow / num_ele)^(1/pow);  % LV = ( 1/n * Simgma_e (g_e)^p )^(1/p)
for ele_idx = 1:num_ele
    LocalVe = LocalVe_All(ele_idx);
    SumLocalVe0 = SumLocalVe0_All(ele_idx);
    Neighbors_idx = Neighbors_ECell{ele_idx};
    for nidx = 1:numel(Neighbors_idx)
        pat_idx = EleInfoMat(Neighbors_idx(nidx),1);
        row_idx = EleInfoMat(Neighbors_idx(nidx),2);
        col_idx = EleInfoMat(Neighbors_idx(nidx),3);
        
        ne_u = num_ele_Mat(pat_idx,1);
        iijj = (col_idx-1)*ne_u + row_idx; % local idx within the patch: pat_idx, for this element
        Ve_All = StiffnessVolMat_PCell{pat_idx,2};
        Weight_All = StiffnessVolMat_PCell{pat_idx,3};
        DNonZeroIdx = DNonZeroIdx_PCell{pat_idx}{iijj};
        NonZeroBasis_GuassPt_PCell = NonZeroBasis_GuassPt_PCell{pat_idx}{iijj};
        
        for kgp = 1:ngst
            DerHDensiijj = DerHDense_All(Neighbors_idx(nidx),kgp);
            NIIJJ = NonZeroBasis_GuassPt_PCell(kgp,:)';
            NodeIdxIIJJ = DNonZeroIdx(kgp,:)';
            DersLocalVol(NodeIdxIIJJ) = DersLocalVol(NodeIdxIIJJ) + LocalVe^(pow-1) / SumLocalVe0 * Ve_All(iijj,kgp) * Weight_All(iijj,kgp) * DerHDensiijj * NIIJJ;
        end
        
    end
end

DersLocalVol = DersLocalVol * ((SumLVe_pow/num_ele)^(1/pow-1)) / num_ele;

end

