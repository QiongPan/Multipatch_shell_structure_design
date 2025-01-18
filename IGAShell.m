/**
 * @file IGAShell.m
 * @author Qiong Pan (PQ2019@mail.ustc.edu.cn)
 * @brief  Code for <<Isogeometric Topology Optimization of Multi-patch Shell Structures>>.
 * Orginal work written by Qiong Pan in Matlab.
 * @date 2024-01-23
 *  Please feel free to contact us with any questions! 
 */
 
function IGAShell()
%%%%%% the main program file that used to optimize the topology of shell structure %%%%%%
% user defined:
%   sType                   - the label of mid-surface
%                           - 1-2-Hyperbolic, 4-torus, 5-table holder, 6-pipes, 3-bird's nest
%   multiPatSurf            - the geometric information of the mid-surface, = NurbsSurface(sType)
%                           - a cell with dim = num_pat * 1
%                           - multiPatSurf(k,:) = {p,q,m,n,U,V,Ctrlpts}
%   num_ele                 - the maximum number of knot spans in a patch edge among all patches for analysis model
%   DctrlElenum             - the maximum number of knot spans in a patch edge among all patches for design model
%   Gload                   - the magnitude of the force
%   loadTypeVec             - the different way to apply force
%   ploadsMat_PCell         - the load applied on the patch
%                           - a matrix with 5 columns
%                           - ploads_Mat(i,:) = (s_i,t_i,fx_i,fy_i,fz_i)
%   fixed_DOF               - fixed control points
%   myFilter                - volume constraint type
%   filter_radius           - a parameter related to local influence radius(Radius) in local volum constraint
%   volfrac                 - the upper bound in volume constraint  
% output(auto saved):
%   xval                    - the value of design variables(coefficients)
%   Design_UV_PCell         - U-,V-knot vectors in design model         (patch related cell)
%   Design_Weights_PCell    - weights in NURBS defined density function (patch related cell)
%   Design_Patch_DOF_PCell  - global index of design variables for each patch
%   num_ele_Mat             - the numbers of knot spans in patch edges for all patches for analysis model
%   interfaceMat            - the interface information for patch edges, a matrix with dim = num_pat * 4
%   Data                    - the iteration information
%   multiPatSurf_           - the geometric information of refined matching patches.
%   save(path, 'xval', 'Design_UV_PCell','Design_Weights_PCell','Design_Patch_DOF_PCell','num_ele_Mat','interfaceMat','Data','multiPatSurf_');

clear all; clc;

%% global constant
myFilter = 'NULL';                      % total volume constraint
% myFilter = 'LocalVolumeConstraint';   % local volume constraint for porous structure design with local influence radius 'Radius':  
                                        % |Pi-Pj| <= Radius = filter_adius * elements_edge_length_average

for kij = 1:5
        global E0 Emin poisson thickness penal beta ng_u ng_v ng_th p q filter_radius;
        E0          = 1;
        poisson     = 0.3;
        thickness   = 1;
        penal       = 05;
        beta        = 2;
        Emin        = 1E-5;
        %       ng_u         = 3;
        %       ng_v         = 3;
        ng_u         = 2;
        ng_v         = 2;
        ng_th        = 3;
        filter_radius = 0;
        
        %% sType = 1; Hyperbolic -- 1 patch
        %         sType = 1;
        %         Gload = 0.1;
        %         Ele_kij = [ 
        %             100 50;
        %             200 50;];
        %         num_ele = Ele_kij(kij,1);
        %         DctrlElenum = Ele_kij(kij,2);
        %% sType = 2; Hyperbolic -- 4 patches
        %         sType = 2;
        %         Gload = 0.1;
        %         Ele_kij = [
        %             64 32;
        %             128 32;];
        %         num_ele = Ele_kij(kij,1);
        %         DctrlElenum = Ele_kij(kij,2);
        %% sType = 4; torus -- 4 patches
        sType = 4;
        Gload = 1;
        Ele_kij = [ 
            175 42;
            175 56;
            175 70;
            175 105;
            175 140;
            ];
        num_ele = Ele_kij(kij,1);
        DctrlElenum = Ele_kij(kij,2);
        %% sType = 5; tablet holder -- 4 patches 
        %         sType = 5;
        %         Gload = 5E-2;
        %         Ele_kij = [ 
        %             50 20;
        %             50 30;];
        %         num_ele = Ele_kij(kij,1);
        %         DctrlElenum = Ele_kij(kij,2);
        %% sType = 6; pipes -- 24 patches
        %         sType = 6;
        %         Gload = 5E-2;
        %         Ele_kij = [ 
        %             50 10;];
        %         num_ele = Ele_kij(kij,1);
        %         DctrlElenum = Ele_kij(kij,2);
        %% sType = 3; bird's nest -- 4 patches
        %         sType = 7;
        %         Gload = 5E-3;
        %         Ele_kij = [
        %             80 40;];
        %         num_ele = Ele_kij(kij,1);
        %         DctrlElenum = Ele_kij(kij,2);
        %%
        if (kij > size(Ele_kij,1))
            break
        end
        multiPatSurf = NurbsSurface(sType); % {p,q,m,n,U,V,Ctrlpts};
        
        %% input a multi-patch surface, extract corner points and select interface candidates based on them
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find interface candidates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Notice: Not applicable to degenerate edges and closed surfaces
        num_pat   = size(multiPatSurf,1); % the number of patches
        CorPtsMat = zeros(num_pat, 4, 3); % corner point np*4*4
        for k = 1:num_pat
            CPs = cell2mat(multiPatSurf(k,7));
            m_k = cell2mat(multiPatSurf(k,3));
            n_k = cell2mat(multiPatSurf(k,4));
            CorPtsMat(k,:,:) = reshape([CPs(1,1,1:3); CPs(m_k+1,1,1:3); CPs(m_k+1,n_k+1,1:3); CPs(1,n_k+1,1:3)],4,3);
        end
        
        pseudo_interface_info = [];	% storage the (patch1,bdr1,patch2,bdr2) information for each interface candidate.
        num_interface_condidate = 0;
        CorPtsMat_ = zeros(num_pat*4,3);
        CorPtsMat_(:,1) = reshape(CorPtsMat(:,:,1)',num_pat*4,1);
        CorPtsMat_(:,2) = reshape(CorPtsMat(:,:,2)',num_pat*4,1);
        CorPtsMat_(:,3) = reshape(CorPtsMat(:,:,3)',num_pat*4,1);
        
        if(num_pat > 1)
            for k = 1:num_pat % patch loop
                for L = 1:4   % bdr loop
                    
                    switch L     % bdr L of patch k
                        case 1
                            CorPt_1 = reshape(CorPtsMat(k,1,:),1,3);
                            CorPt_2 = reshape(CorPtsMat(k,2,:),1,3);
                        case 2
                            CorPt_1 = reshape(CorPtsMat(k,2,:),1,3);
                            CorPt_2 = reshape(CorPtsMat(k,3,:),1,3);
                        case 3
                            CorPt_1 = reshape(CorPtsMat(k,4,:),1,3);
                            CorPt_2 = reshape(CorPtsMat(k,3,:),1,3);
                        case 4
                            CorPt_1 = reshape(CorPtsMat(k,1,:),1,3);
                            CorPt_2 = reshape(CorPtsMat(k,4,:),1,3);
                    end
                    
                    diffMat_1 = repmat(CorPt_1,num_pat*4,1) - CorPtsMat_;
                    distMat_1 = sum(diffMat_1.*diffMat_1, 2);
                    idx_CorPts = find(distMat_1 == 0);
                    idx_CorPts = setdiff(idx_CorPts,reshape((0:k-1)*4+(1:4)',[],1));    % other patches(haven't been checked) that shares the same point CorPt_1.
                    
                    for a = 1:length(idx_CorPts)
                        idx = idx_CorPts(a);
                        idx_pat =ceil(idx / 4); % patch k_
                        idx1 = mod(idx, 4);     % point idx in patch k_
                        
                        if(idx1 == 0)
                            idx1 = 4;
                        end
                        
                        switch idx1   % find another point in patch k_, i.e. find the bdr L_ in patch k_.
                            case 1
                                CorPt_ = reshape(CorPtsMat(idx_pat,2,:),1,3);
                                if(norm(CorPt_-CorPt_2) == 0)   % L_ = 1
                                    num_interface_condidate = num_interface_condidate + 1;
                                    pseudo_interface_info = [pseudo_interface_info; [k, L, idx_pat, 1]];
                                else                            % L_ = 4
                                    CorPt_ = reshape(CorPtsMat(idx_pat,4,:),1,3);
                                    if(norm(CorPt_-CorPt_2) == 0)
                                        num_interface_condidate = num_interface_condidate + 1;
                                        pseudo_interface_info = [pseudo_interface_info; [k, L, idx_pat, 4]];
                                    end
                                end
                            case 2
                                CorPt_ = reshape(CorPtsMat(idx_pat,3,:),1,3);
                                if(norm(CorPt_-CorPt_2) == 0)   % L_ = 2
                                    num_interface_condidate = num_interface_condidate + 1;
                                    pseudo_interface_info = [pseudo_interface_info; [k, L, idx_pat, 2]];
                                else                            % L_ = 1
                                    CorPt_ = reshape(CorPtsMat(idx_pat,1,:),1,3);
                                    if(norm(CorPt_-CorPt_2) == 0)
                                        num_interface_condidate = num_interface_condidate + 1;
                                        pseudo_interface_info = [pseudo_interface_info; [k, L, idx_pat, 1]];
                                    end
                                end
                            case 3
                                CorPt_ = reshape(CorPtsMat(idx_pat,4,:),1,3);
                                if(norm(CorPt_-CorPt_2) == 0)   % L_ = 3
                                    num_interface_condidate = num_interface_condidate + 1;
                                    pseudo_interface_info = [pseudo_interface_info; [k, L, idx_pat, 3]];
                                else                            % L_ = 2
                                    CorPt_ = reshape(CorPtsMat(idx_pat,2,:),1,3);
                                    if(norm(CorPt_-CorPt_2) == 0)
                                        num_interface_condidate = num_interface_condidate + 1;
                                        pseudo_interface_info = [pseudo_interface_info; [k, L, idx_pat, 2]];
                                    end
                                    
                                end
                            case 4
                                CorPt_ = reshape(CorPtsMat(idx_pat,3,:),1,3);
                                if(norm(CorPt_-CorPt_2) == 0)   % L_ = 3
                                    num_interface_condidate = num_interface_condidate + 1;
                                    pseudo_interface_info = [pseudo_interface_info; [k, L, idx_pat, 3]];
                                else                            % L_ = 4
                                    CorPt_ = reshape(CorPtsMat(idx_pat,1,:),1,3);
                                    if(norm(CorPt_-CorPt_2) == 0)
                                        num_interface_condidate = num_interface_condidate + 1;
                                        pseudo_interface_info = [pseudo_interface_info; [k, L, idx_pat, 4]];
                                    end
                                end
                        end
                    end
                    
                end
            end
        end
        
        % Traversing the same corner in this way may find 4 or more patches sharing this point, leading to these patches be traversed at least twice.
        % For example: L1 and L4 both contain P1, and the patches adjacent to L1 and L4 will be traversed twice.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% refine candidate adjacent patches for further determination of interfaces
        %%%%%%%%%%%%%%% unify the representation of the curve by refinement: order elevation, knots insertion %%%%%%%%%%%%%%%%
        % Notice：not suitable for adjacent patches with opposite parameter directions along the interface
        % The representation corresponding to the same NURBS curve if the degree, knot vectors, and control points are the same,
        % and the weights of the two sets of basis functions differ by a constant factor.
        bool_interface = ones(num_interface_condidate,1);
        multiPatSurf_ = multiPatSurf;
        if(num_pat > 1)
            a = 1;
            while (a <= num_interface_condidate)
                if(bool_interface(a)==0)
                    a = a + 1;
                    continue;
                end
                
                idx_pat1 = pseudo_interface_info(a,1); idx_edge1 = pseudo_interface_info(a,2);
                idx_pat2 = pseudo_interface_info(a,3); idx_edge2 = pseudo_interface_info(a,4);
                
                p1 = cell2mat(multiPatSurf_(idx_pat1,1));
                q1 = cell2mat(multiPatSurf_(idx_pat1,2));
                m1 = cell2mat(multiPatSurf_(idx_pat1,3));
                n1 = cell2mat(multiPatSurf_(idx_pat1,4));
                U1 = cell2mat(multiPatSurf_(idx_pat1,5));
                V1 = cell2mat(multiPatSurf_(idx_pat1,6));
                intfPat1_CtrlPts_ = cell2mat(multiPatSurf_(idx_pat1,7));
                switch idx_edge1
                    case 1
                        intf_deg_1 = p1;
                        intf_knots_1 = U1;
                        intf_CtrlPts_1 = reshape(intfPat1_CtrlPts_(1:m1+1,1,:),m1+1,4);
                    case 2
                        intf_deg_1 = q1;
                        intf_knots_1 = V1;
                        intf_CtrlPts_1 = reshape(intfPat1_CtrlPts_(m1+1,1:n1+1,:),n1+1,4);
                    case 3
                        intf_deg_1 = p1;
                        intf_knots_1 = U1;
                        intf_CtrlPts_1 = reshape(intfPat1_CtrlPts_(1:m1+1,n1+1,:),m1+1,4);
                    case 4
                        intf_deg_1 = q1;
                        intf_knots_1 = V1;
                        intf_CtrlPts_1 = reshape(intfPat1_CtrlPts_(1,1:n1+1,:),n1+1,4);
                end
                
                p2 = cell2mat(multiPatSurf_(idx_pat2,1));
                q2 = cell2mat(multiPatSurf_(idx_pat2,2));
                m2 = cell2mat(multiPatSurf_(idx_pat2,3));
                n2 = cell2mat(multiPatSurf_(idx_pat2,4));
                U2 = cell2mat(multiPatSurf_(idx_pat2,5));
                V2 = cell2mat(multiPatSurf_(idx_pat2,6));
                intfPat2_CtrlPts_ = cell2mat(multiPatSurf_(idx_pat2,7));
                switch idx_edge2
                    case 1
                        intf_deg_2 = p2;
                        intf_knots_2 = U2;
                        intf_CtrlPts_2 = reshape(intfPat2_CtrlPts_(1:m2+1,1,:),m2+1,4);
                    case 2
                        intf_deg_2 = q2;
                        intf_knots_2 = V2;
                        intf_CtrlPts_2 = reshape(intfPat2_CtrlPts_(m2+1,1:n2+1,:),n2+1,4);
                    case 3
                        intf_deg_2 = p2;
                        intf_knots_2 = U2;
                        intf_CtrlPts_2 = reshape(intfPat2_CtrlPts_(1:m2+1,n2+1,:),m2+1,4);
                    case 4
                        intf_deg_2 = q2;
                        intf_knots_2 = V2;
                        intf_CtrlPts_2 = reshape(intfPat2_CtrlPts_(1,1:n2+1,:),n2+1,4);
                end
                
                if(intf_deg_1 == intf_deg_2) % with the same degree
                    if(length(intf_knots_1)==length(intf_knots_2) && sum(intf_knots_1~=intf_knots_2)==0) % with the same knot vector
                        if(sum(abs(intf_CtrlPts_1(:,1:3)-intf_CtrlPts_2(:,1:3)),'all') < 1E-2 && isequal(intf_CtrlPts_1(:,4)/norm(intf_CtrlPts_1(:,4)),intf_CtrlPts_2(:,4)/norm(intf_CtrlPts_2(:,4)))) % with the same control points
                            [rows_1, ~] = find(pseudo_interface_info(:,1) == idx_pat1 & pseudo_interface_info(:,2) == idx_edge1);
                            [rows_2, ~] = find(pseudo_interface_info(:,1) == idx_pat2 & pseudo_interface_info(:,2) == idx_edge2);
                            [rows_3, ~] = find(pseudo_interface_info(:,3) == idx_pat1 & pseudo_interface_info(:,4) == idx_edge1);
                            [rows_4, ~] = find(pseudo_interface_info(:,3) == idx_pat2 & pseudo_interface_info(:,4) == idx_edge2);
                            rows = setdiff(union(union(rows_1,rows_2),union(rows_3, rows_4)),a);
                            bool_interface(rows) = 0;
                            a = a + 1; % go to next patch
                            multiPatSurf_(idx_pat1,:) = multiPatSurf(idx_pat1,:);
                            multiPatSurf_(idx_pat2,:) = multiPatSurf(idx_pat2,:);
                        else % with the same degree and knot vector but with different control points
                            bool_interface(a) = 0;
                        end
                    else % with different degree
                        % merge knot vector and insertknots mutually
                        UniqueKnots = union(intf_knots_1,intf_knots_2);
                        KnotsInsert_1 = [];
                        KnotsInsert_2 = [];
                        for b = 1:length(UniqueKnots)
                            uv = UniqueKnots(b);
                            mult_1 = FindMultiplicity(uv, intf_knots_1);
                            mult_2 = FindMultiplicity(uv, intf_knots_2);
                            mult_diff = mult_1 - mult_2;
                            if(mult_diff > 0)
                                KnotsInsert_2 = [KnotsInsert_2; ones(mult_diff,1)*uv];
                            else
                                KnotsInsert_1 = [KnotsInsert_1; ones(-mult_diff,1)*uv];
                            end
                        end
                        
                        intfPat1_HCtrlPts_ = intfPat1_CtrlPts_;
                        intfPat1_HCtrlPts_(:,:,1:3) = intfPat1_HCtrlPts_(:,:,1:3).*intfPat1_HCtrlPts_(:,:,4);
                        if(mod(idx_edge1,2) == 1) %if(idx_edge1 == 1 || idx_edge1 == 3)   %Insert knots in u-direction
                            [U1, V1, intfPat1_HCtrlPts_] = KnotInsert(KnotsInsert_1,[],U1,V1,intfPat1_HCtrlPts_);
                        else                                                               %Insert knots in v-direction
                            [U1, V1, intfPat1_HCtrlPts_] = KnotInsert([],KnotsInsert_1,U1,V1,intfPat1_HCtrlPts_);
                        end
                        
                        intfPat2_HCtrlPts_ = intfPat2_CtrlPts_;
                        intfPat2_HCtrlPts_(:,:,1:3) = intfPat2_HCtrlPts_(:,:,1:3).*intfPat2_HCtrlPts_(:,:,4);
                        if(mod(idx_edge2,2) == 1) %if(idx_edge2 == 1 || idx_edge2 == 3)
                            [U2, V2, intfPat2_HCtrlPts_] = KnotInsert(KnotsInsert_2,[],U2,V2,intfPat2_HCtrlPts_);
                        else
                            [U2, V2, intfPat2_HCtrlPts_] = KnotInsert([],KnotsInsert_2,U2,V2,intfPat2_HCtrlPts_);
                        end
                        
                        [m1,n1,~] = size(intfPat1_HCtrlPts_);
                        [m2,n2,~] = size(intfPat2_HCtrlPts_);
                        intfPat1_CtrlPts_ = intfPat1_HCtrlPts_;
                        intfPat1_CtrlPts_(:,:,1:3) = intfPat1_HCtrlPts_(:,:,1:3)./intfPat1_HCtrlPts_(:,:,4);
                        intfPat2_CtrlPts_ = intfPat2_HCtrlPts_;
                        intfPat2_CtrlPts_(:,:,1:3) = intfPat2_HCtrlPts_(:,:,1:3)./intfPat2_HCtrlPts_(:,:,4);
                        multiPatSurf_(idx_pat1,:) = {p1,q1,m1-1,n1-1,U1,V1,intfPat1_CtrlPts_};
                        multiPatSurf_(idx_pat2,:) = {p2,q2,m2-1,n2-1,U2,V2,intfPat2_CtrlPts_};
                    end
                elseif(intf_deg_1 > intf_deg_2)
                    intfPat2_HCtrlPts_ = intfPat2_CtrlPts_;
                    intfPat2_HCtrlPts_(:,:,1:3) = intfPat2_HCtrlPts_(:,:,1:3).*intfPat2_HCtrlPts_(:,:,4);
                    if(mod(idx_edge2,2) == 1) %if(idx_edge2 == 1 || idx_edge2 == 3)  % Direction = 'U-Direction';
                        [intfPat2_HCtrlPts_,U2,V2] = DegreeElevateSurface(p2,q2,U2,V2,intfPat2_HCtrlPts_,intf_deg_1 - intf_deg_2,0);
                        p2 = intf_deg_1;
                    else
                        [intfPat2_HCtrlPts_,U2,V2] = DegreeElevateSurface(p2,q2,U2,V2,intfPat2_HCtrlPts_,0,intf_deg_1 - intf_deg_2);
                        q2 = intf_deg_1;
                    end
                    [m2,n2,~] = size(intfPat2_HCtrlPts_);
                    intfPat2_CtrlPts_ = intfPat2_HCtrlPts_;
                    intfPat2_CtrlPts_(:,:,1:3) = intfPat2_HCtrlPts_(:,:,1:3)./intfPat2_HCtrlPts_(:,:,4);
                    multiPatSurf_(idx_pat2,:) = {p2,q2,m2-1,n2-1,U2,V2,intfPat2_CtrlPts_};
                else % (intf_deg_1 < intf_deg_2)
                    intfPat1_HCtrlPts_ = intfPat1_CtrlPts_;
                    intfPat1_HCtrlPts_(:,:,1:3) = intfPat1_HCtrlPts_(:,:,1:3).*intfPat1_HCtrlPts_(:,:,4);
                    if(mod(idx_edge1,2) == 1) %if(idx_edge1 == 1 || idx_edge1 == 3)
                        [intfPat1_HCtrlPts_,U1,V1] = DegreeElevateSurface(p1,q1,U1,V1,intfPat1_HCtrlPts_,intf_deg_2 - intf_deg_1,0);
                        p1 = intf_deg_2;
                    else
                        [intfPat1_HCtrlPts_,U1,V1] = DegreeElevateSurface(p1,q1,U1,V1,intfPat1_HCtrlPts_,intf_deg_2 - intf_deg_1,0);
                        q1 = intf_deg_2;
                    end
                    [m1,n1,~] = size(intfPat1_HCtrlPts_);
                    intfPat1_CtrlPts_ = intfPat1_HCtrlPts_;
                    intfPat1_CtrlPts_(:,:,1:3) = intfPat1_HCtrlPts_(:,:,1:3)./intfPat1_HCtrlPts_(:,:,4);
                    multiPatSurf_(idx_pat1,:) = {p1,q1,m1-1,n1-1,U1,V1,intfPat1_CtrlPts_};
                end
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% distinguish true interface
        num_interface = sum(bool_interface);
        interfaceMat = zeros(num_pat, 4); % interface indicies for four edges
        interface_idx_in_pseudoInfo = find(bool_interface);
        for int = 1:num_interface
            idx = interface_idx_in_pseudoInfo(int);
            patch_1 = pseudo_interface_info(idx,1); edge_1 = pseudo_interface_info(idx,2);
            patch_2 = pseudo_interface_info(idx,3); edge_2 = pseudo_interface_info(idx,4);
            interfaceMat(patch_1,edge_1) = int;
            interfaceMat(patch_2,edge_2) = int;
        end
        
        %% compute length of edges: \int_0^1 \sqrt{x'^2 + y'^2 + z'^2} dt
        pat_len_wid_Mat = zeros(num_pat,4);
        for k = 1:num_pat
            p_k = cell2mat(multiPatSurf(k,1));
            q_k = cell2mat(multiPatSurf(k,2));
            U_k = cell2mat(multiPatSurf(k,5));
            V_k = cell2mat(multiPatSurf(k,6));
            Ctrlpts_k = cell2mat(multiPatSurf(k,7));
            HCtrlpts_k = Ctrlpts_k;
            HCtrlpts_k(:,:,1:3) = HCtrlpts_k(:,:,1:3).*HCtrlpts_k(:,:,4);
            
            [c1,c2,c3,c4] = ComputeCircumference(p_k,q_k,U_k,V_k,HCtrlpts_k);
            pat_len_wid_Mat(k,:) = [c1,c2,c3,c4];
        end
        
        %% find related edges/interfaces and unify labeling of edges for subsequent global refinement
        %%%%%%%%%%%%%%%%%%%%%%%%%%% unify labeling of related edges for subsequent global refinement %%%%%%%%%%%%%%%%%%%%%%%%%
        interfaceMat_ = interfaceMat;
        for int = 1:num_interface
            [rows, cols] = find(interfaceMat_ == int);
            if(numel(rows) == 0)
                continue
            else
                while(~isempty(rows))
                    ii = 1;
                    row = rows(ii);
                    col = cols(ii);
                    col_temp = mod(col+2,4);
                    if(col_temp == 0)
                        col_temp = 4;
                    end
                    int_temp = interfaceMat_(row,col_temp);
                    if(int_temp~=0 && int_temp~=int) % the opposite side of the adjacent patch
                        [rows_, cols_] = find(interfaceMat_ == int_temp);
                        [~, ridx] = setdiff(rows_,row);
                        row_ = rows_(ridx);
                        col_ = cols_(ridx);
                        rows(ii) = row_;
                        cols(ii) = col_;
                        interfaceMat_(row_,col_) = int;
                        interfaceMat_(row,col_temp) = int;
                    else
                        rows = rows(ii+1:end);
                        cols = cols(ii+1:end);
                    end
                end
                
            end
        end
        
        % unify labeling of all sides: some patches only have adjacent patches on one side but not on the opposite side
        [rows, cols] = find(interfaceMat_ == 0);
        for ii = 1:numel(rows)
            row = rows(ii);
            col = cols(ii);
            col_temp = mod(col+2,4);
            if(col_temp == 0)
                col_temp = 4;
            end
            interfaceMat_(row,col) = interfaceMat_(row,col_temp);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% unify the degree along the interface for all adjacent patches
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pat_deg_Mat = zeros(num_pat,2); % num_pat * ([num_uknots * 2], [num_vknots * 2])
        for k = 1:num_pat
            p_k = cell2mat(multiPatSurf(k,1));
            q_k = cell2mat(multiPatSurf(k,2));
            pat_deg_Mat(k,:) = [p_k, q_k];
        end
        
        % find union degree of interface for comformming patches
        UnionPatDeg_Mat = zeros(num_interface,1);
        for int = 1:num_interface
            [rows, cols] = find(interfaceMat_ == int);
            if(numel(rows) == 0)
                continue
            else
                row = rows(1);
                col = cols(1);
                if(mod(col,2) == 1)
                    Deg_1 = pat_deg_Mat(row,1);
                else
                    Deg_1 = pat_deg_Mat(row,2);
                end
                for ii = 2: numel(rows)
                    row = rows(ii);
                    col = cols(ii);
                    if(mod(col,2) == 1)
                        Deg_2 = pat_deg_Mat(row,1);
                    else
                        Deg_2 = pat_deg_Mat(row,2);
                    end
                    Deg_1 = max(Deg_1, Deg_2);
                end
                
                UnionPatDeg_Mat(int) = Deg_1;
            end
        end
        
        % elevate degree for patches accross interface
        multiPatSurf_ = multiPatSurf;
        for int = 1:num_interface
            [rows, cols] = find(interfaceMat_ == int);
            if(numel(rows) == 0)
                continue;
            else
                Deg = UnionPatDeg_Mat(int);
                
                for ii = 1: numel(rows)
                    row = rows(ii);
                    col = cols(ii);
                    
                    p_ = cell2mat(multiPatSurf_(row,1));
                    q_ = cell2mat(multiPatSurf_(row,2));
                    U_ = cell2mat(multiPatSurf_(row,5));
                    V_ = cell2mat(multiPatSurf_(row,6));
                    Pat_CtrlPts_ = cell2mat(multiPatSurf_(row,7));
                    
                    if(mod(col,2) == 1)
                        Deg_ = p_;
                    else
                        Deg_ = q_;
                    end
                    
                    Pat_HCtrlPts_ = Pat_CtrlPts_;
                    Pat_HCtrlPts_(:,:,1:3) = Pat_HCtrlPts_(:,:,1:3).*Pat_HCtrlPts_(:,:,4);
                    if(mod(col,2) == 1)
                        [Pat_HCtrlPts_, U_, V_] = DegreeElevateSurface(p_,q_,U_,V_,Pat_HCtrlPts_,Deg - Deg_,0);
                        p_ = Deg;
                    else
                        [Pat_HCtrlPts_, U_, V_] = DegreeElevateSurface(p_,q_,U_,V_,Pat_HCtrlPts_,0,Deg - Deg_);
                        q_ = Deg;
                    end
                    
                    [m_,n_,~] = size(Pat_HCtrlPts_);
                    Pat_CtrlPts_ = Pat_HCtrlPts_;
                    Pat_CtrlPts_(:,:,1:3) = Pat_HCtrlPts_(:,:,1:3)./Pat_HCtrlPts_(:,:,4);
                    multiPatSurf_(row,:) = {p_,q_,m_-1,n_-1,U_,V_,Pat_CtrlPts_};
                end
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% unify the knot vector along the interface for all adjacent patches
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pat_knot_mult_Cell = cell(num_pat,2); % num_pat * ([num_uknots * 2], [num_vknots * 2])
        for k = 1:num_pat
            U_k = cell2mat(multiPatSurf_(k,5));
            V_k = cell2mat(multiPatSurf_(k,6));
            Uk_diffknots = reshape(unique(U_k),[],1);
            Vk_diffknots = reshape(unique(V_k),[],1);
            Uk_mult = zeros(numel(Uk_diffknots),1);
            Vk_mult = zeros(numel(Vk_diffknots),1);
            for i = 1:numel(Uk_diffknots)
                ui = Uk_diffknots(i);
                Uk_mult(i) = FindMultiplicity(ui, U_k);
            end
            for j = 1:numel(Vk_diffknots)
                vj = Vk_diffknots(j);
                Vk_mult(j) = FindMultiplicity(vj, V_k);
            end
            pat_knot_mult_Cell(k,:) = {[Uk_diffknots,Uk_mult],[Vk_diffknots,Vk_mult]};
        end
        
        % find union knots and multiplicities for each interface group.
        UnionPatKM_Cell = cell(num_interface,1); % num_pat * ([num_knots * 2])
        for int = 1:num_interface
            [rows, cols] = find(interfaceMat_ == int);
            if(numel(rows) == 0)
                UnionPatKM_Cell(int) = {[]};
                continue
            else
                row = rows(1);
                col = cols(1);
                if(mod(col,2) == 1)
                    Knot_Mult_1 = cell2mat(pat_knot_mult_Cell(row,1));
                else
                    Knot_Mult_1 = cell2mat(pat_knot_mult_Cell(row,2));
                end
                for ii = 2: numel(rows)
                    row = rows(ii);
                    col = cols(ii);
                    if(mod(col,2) == 1)
                        Knot_Mult_2 = cell2mat(pat_knot_mult_Cell(row,1));
                    else
                        Knot_Mult_2 = cell2mat(pat_knot_mult_Cell(row,2));
                    end
                    
                    Knots = union(Knot_Mult_1(:,1),Knot_Mult_2(:,1));
                    Mults = zeros(numel(Knots),1);
                    for jj = 1:numel(Knots)
                        uv_jj = Knots(jj);
                        idx_1 = (Knot_Mult_1(:,1) == uv_jj);
                        idx_2 = (Knot_Mult_2(:,1) == uv_jj);
                        mult_1 = Knot_Mult_1(idx_1, 2);
                        mult_2 = Knot_Mult_2(idx_2, 2);
                        if (isempty(mult_1))
                            mult_1 = 0;
                        end
                        if (isempty(mult_2))
                            mult_2 = 0;
                        end
                        
                        Mults(jj) = max(mult_1, mult_2);
                    end
                    Knot_Mult_1 = [Knots, Mults];
                end
                
                UnionPatKM_Cell(int) = {Knot_Mult_1};
            end
        end
        
        % insert knots for patches to achive uniform Knots_Multiplicities in a interface group.
        for int = 1:num_interface
            [rows, cols] = find(interfaceMat_ == int);
            if(numel(rows) == 0)
                continue;
            else
                Knot_Mult = cell2mat(UnionPatKM_Cell(int));
                for ii = 1: numel(rows)
                    row = rows(ii);
                    col = cols(ii);
                    
                    p_ = cell2mat(multiPatSurf_(row,1));
                    q_ = cell2mat(multiPatSurf_(row,2));
                    U_ = cell2mat(multiPatSurf_(row,5));
                    V_ = cell2mat(multiPatSurf_(row,6));
                    Pat_CtrlPts_ = cell2mat(multiPatSurf_(row,7));
                    
                    KnotsInsert_ = [];
                    if(mod(col,2) == 1)
                        Knot_Mult_ = cell2mat(pat_knot_mult_Cell(row,1));
                    else
                        Knot_Mult_ = cell2mat(pat_knot_mult_Cell(row,2));
                    end
                    for jj = 1:size(Knot_Mult,1)
                        uv_jj = Knot_Mult(jj, 1);
                        mult = Knot_Mult(jj,2);
                        idx_ = find(Knot_Mult_(:,1) == uv_jj);
                        if(isempty(idx_))
                            mult_ = 0;
                        else
                            mult_ = Knot_Mult_(idx_,2);
                        end
                        mult_diff = mult - mult_;
                        KnotsInsert_ = [KnotsInsert_; ones(mult_diff,1)*uv_jj];
                    end
                    
                    Pat_HCtrlPts_ = Pat_CtrlPts_;
                    Pat_HCtrlPts_(:,:,1:3) = Pat_HCtrlPts_(:,:,1:3).*Pat_HCtrlPts_(:,:,4);
                    if(mod(col,2) == 1)
                        [U_, V_, Pat_HCtrlPts_] = KnotInsert(KnotsInsert_,[],U_,V_,Pat_HCtrlPts_);
                        pat_knot_mult_Cell(row,1) = {Knot_Mult};
                    else
                        [U_, V_, Pat_HCtrlPts_] = KnotInsert([],KnotsInsert_,U_,V_,Pat_HCtrlPts_);
                        pat_knot_mult_Cell(row,2) = {Knot_Mult};
                    end
                    
                    [m_,n_,~] = size(Pat_HCtrlPts_);
                    Pat_CtrlPts_ = Pat_HCtrlPts_;
                    Pat_CtrlPts_(:,:,1:3) = Pat_HCtrlPts_(:,:,1:3)./Pat_HCtrlPts_(:,:,4);
                    multiPatSurf_(row,:) = {p_,q_,m_-1,n_-1,U_,V_,Pat_CtrlPts_};
                end
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% subdivide elements according to length of edges
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_ele_Mat = zeros(num_pat, 2);
        Mall = max(pat_len_wid_Mat,[],'all');
        for k = 1: num_pat
            len = max(pat_len_wid_Mat(k,1),pat_len_wid_Mat(k,3));
            wid = max(pat_len_wid_Mat(k,2),pat_len_wid_Mat(k,4));
            num_ele_Mat(k,1) = round(len/Mall*num_ele);
            num_ele_Mat(k,2) = round(wid/Mall*num_ele);
        end
        
        % uniform element subdivision for interfaces
        num_ele_Mat_ = repmat(num_ele_Mat,1,2);
        for int = 1:num_interface
            [rows, cols] = find(interfaceMat_ == int);
            if(numel(rows) == 0)
                continue
            else
                ind = sub2ind([num_pat,4],rows,cols);
                max_ele_num = max(num_ele_Mat_(ind));
                for ii = 1: numel(rows)
                    row = rows(ii);
                    col = cols(ii);
                    if(mod(col,2) == 1)
                        num_ele_Mat(row,1) = max_ele_num;
                    else
                        num_ele_Mat(row,2) = max_ele_num;
                    end
                end
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Compute stiffness matrix and volume for patches on Gauss integration points: StiffnessVolWeight_PCell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparation for Gauss integration(computed only once) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        UEle_Cell = cell(num_pat,1);
        VEle_Cell = cell(num_pat,1);
        ANonZeroIdx_PCell        = cell(num_pat,1);
        multiPatSurf_Analysis   = multiPatSurf_;
        StiffnessVolWeight_PCell    = cell(num_pat,3); %{Ke_GuassPt_ECell, Ve_GuassPt_All, Weight_GuassPt_All}
        
        for k = 1:num_pat
            p = cell2mat(multiPatSurf_(k,1));
            q = cell2mat(multiPatSurf_(k,2));
            U = cell2mat(multiPatSurf_(k,5));
            V = cell2mat(multiPatSurf_(k,6));
            Ctrlpts = cell2mat(multiPatSurf_(k,7));
            HCtrlpts = Ctrlpts;
            HCtrlpts(:,:,1:3) = HCtrlpts(:,:,1:3).*HCtrlpts(:,:,4);
            
            ne_s_more = num_ele_Mat(k,1) - (numel(unique(U))-1);
            ne_t_more = num_ele_Mat(k,2) - (numel(unique(V))-1);
            
            switch sType
                case 4
                    tmp = floor(ne_s_more/2);
                    UEle = [linspace(0,0.5,tmp+2) linspace(0.5,1,tmp+2)];
                    tmp = floor(ne_t_more/2);
                    VEle = [linspace(0,0.5,tmp+2) linspace(0.5,1,tmp+2)];
                otherwise
                    UEle = linspace(U(1),U(end),ne_s_more+2);
                    VEle = linspace(V(1),V(end),ne_t_more+2);
            end
            
            uins = setdiff(UEle, U);
            vins = setdiff(VEle, V);
            
            [UAnalysis, VAnalysis, HCtrlpts] = KnotInsert(uins,vins,U,V,HCtrlpts);
            Ctrlpts             = HCtrlpts;
            Ctrlpts(:,:,1:3)  	= Ctrlpts(:,:,1:3)./Ctrlpts(:,:,4);
            [mAnalysis,nAnalysis,~]    = size(Ctrlpts);
            mAnalysis                  = mAnalysis - 1;
            nAnalysis                  = nAnalysis - 1;
            
            UEle = unique(UAnalysis);
            VEle = unique(VAnalysis);
            ne_u = length(UEle) - 1;
            ne_v = length(VEle) - 1;
            num_ele_Mat(k,:) = [ne_u, ne_v];
            
            [Ke_GuassPt_ECell, Ve_GuassPt_All,Weight_GuassPt_All] = ComputeKV(mAnalysis,nAnalysis,p,q,UAnalysis,VAnalysis,UEle,VEle,Ctrlpts,ne_u,ne_v);
            multiPatSurf_Analysis(k,:) = {p,q,mAnalysis, nAnalysis,UAnalysis,VAnalysis,Ctrlpts};
            StiffnessVolWeight_PCell(k,:) = {Ke_GuassPt_ECell, Ve_GuassPt_All, Weight_GuassPt_All};
            
            % element node (in analysis model)
            UEle_Cell(k) = {UEle};
            VEle_Cell(k) = {VEle};
            
            % indices of non-zero basis (in analysis model) (ne_u*ne_v,(p+1)*(q+1));
            ANonZeroIdx = ComputeNonZeroBasis(p,q,mAnalysis,nAnalysis,UAnalysis,VAnalysis,UEle,VEle,ne_u,ne_v); % dim = (ne_u*ne_v,(p+1)*(q+1));
            ANonZeroIdx_PCell(k) = {ANonZeroIdx};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Compute force vector F
        %%%%%%%%%%%%% The force on the interface is equally divided into two patches and calculated separately %%%%%%%%%%%%%%%
        loadTypeVec = zeros(num_pat,1);
        ploadsMat_PCell = cell(num_pat,1);
        switch sType
            case 1 	% hyperbolic 1 patch
                loadTypeVec = 3;
                ploadsMat_PCell(1) = {[0.5 0.5 0 0 -Gload]};
            case 2 	% hyperbolic 4 patches
                loadTypeVec = [2 3 1 1]';
                ploadsMat_PCell(1) = {[0          0.5/0.6  0  0 -Gload;]};
                ploadsMat_PCell(2) = {[0.15/0.65  0.5/0.6  0  0 -Gload;]};
                ploadsMat_PCell(3) = {[0.15/0.65  0        0  0 -Gload;]};
                ploadsMat_PCell(4) = {[]};
            case 4 	% torus
                fload = Gload/2/sqrt(2);
                ploadsMat_PCell(1) = {[ 0.25 0   -fload fload 0;
                    0.75 0   -fload -fload 0]};
                ploadsMat_PCell(2) = {[ 0.25 1   -fload fload 0;
                    0.75 1  -fload -fload 0]};
                ploadsMat_PCell(3) = {[ 0.25   1   fload -fload 0;
                    0.75 1   fload fload 0]};
                ploadsMat_PCell(4) = {[ 0.25   0   fload -fload 0;
                    0.75 0   fload fload 0]};
                %             case 4
                %                 ploadsMat_PCell(1) = {[0   0   0 Gload/4 0;
                %                                      0.5 0  -Gload/2 0 0;
                %                                      1   0  0 -Gload/4 0]};
                %                 ploadsMat_PCell(2) = {[0   1   0 Gload/4 0;
                %                                      0.5 1  -Gload/2 0 0;
                %                                      1   1  0 -Gload/4 0]};
                %                 ploadsMat_PCell(3) = {[0   1   0 -Gload/4 0;
                %                                      0.5 1   Gload/2  0 0;
                %                                      1   1   0  Gload/4 0]};
                %                 ploadsMat_PCell(4) = {[0   0   0 -Gload/4 0;
                %                                      0.5 0   Gload/2 0 0;
                %                                      1   0  0 Gload/4 0]};
            case 5  % table holder
                ploadsMat_PCell(1) = {[ 0 0 0 0 -Gload/4;
                    1 1 0 0 -Gload]};
                ploadsMat_PCell(2) = {[ 0 0 0 0 -Gload/4;
                    1 1 0 0 -Gload]};
                ploadsMat_PCell(3) = {[ 0 0 0 0 -Gload/4;
                    1 1 0 0 -Gload]};
                ploadsMat_PCell(4) = {[ 0 0 0 0 -Gload/4;
                    1 1 0 0 -Gload]};
            case 6  % pipes
                % computed by 3D rotation
                loadTypeVec(:) = 1;
                num_forcelines = 3;
                theta = linspace(0, pi/2, num_forcelines)';
                alpha_x_axis = pi*[3/2; 1/6; 5/6];
                for pipe_idx = 1:3
                    alpha = alpha_x_axis(pipe_idx);
                    loadline_s = linspace(0, 1, num_forcelines)';
                    forceX =  Gload * sin(theta) * sin(alpha);
                    forceY = -Gload * sin(theta) * cos(alpha);
                    forceZ = -Gload * cos(theta);
                    prevoius_patch_idx = (pipe_idx-1)*8;
                    ploadsMat_PCell(prevoius_patch_idx+1) = {[loadline_s  zeros(num_forcelines,1)      forceX   forceY   forceZ]};
                    ploadsMat_PCell(prevoius_patch_idx+2) = {[loadline_s  zeros(num_forcelines,1)      forceX   forceY  -forceZ]};
                    ploadsMat_PCell(prevoius_patch_idx+3) = {[loadline_s  zeros(num_forcelines,1)     -forceX  -forceY  -forceZ]};
                    ploadsMat_PCell(prevoius_patch_idx+4) = {[loadline_s  zeros(num_forcelines,1)     -forceX  -forceY   forceZ]};
                    
                    ploadsMat_PCell(prevoius_patch_idx+5) = {[loadline_s  zeros(num_forcelines,1)      forceX   forceY   forceZ]};
                    ploadsMat_PCell(prevoius_patch_idx+6) = {[loadline_s  zeros(num_forcelines,1)      forceX   forceY  -forceZ]};
                    ploadsMat_PCell(prevoius_patch_idx+7) = {[loadline_s  zeros(num_forcelines,1)     -forceX  -forceY  -forceZ]};
                    ploadsMat_PCell(prevoius_patch_idx+8) = {[loadline_s  zeros(num_forcelines,1)     -forceX  -forceY   forceZ]};
                end
                for k = 1:num_pat
                    ploadVec_Mat = ploadsMat_PCell{k};
                    ploadVec_Mat([1,num_forcelines],3:5) = ploadVec_Mat([1,num_forcelines],3:5) / 2;
                    ploadsMat_PCell(k) = {ploadVec_Mat};
                end
                for k = 2:4:num_pat
                    ploadVec_Mat = ploadsMat_PCell{k};
                    ploadVec_Mat(1,:) = [];
                    ploadsMat_PCell(k) = {ploadVec_Mat};
                end
                for k = 3:4:num_pat
                    ploadVec_Mat = ploadsMat_PCell{k};
                    ploadVec_Mat(1,:) = [];
                    ploadsMat_PCell(k) = {ploadVec_Mat};
                end
            case 3
                loadTypeVec = ones(1,num_pat);
                ploadsMat_PCell(1) = {[1 0 0 0 -Gload]};
                ploadsMat_PCell(2) = {[1 0 0 0 -Gload]};
                ploadsMat_PCell(3) = {[1 0 0 0 -Gload]};
                ploadsMat_PCell(4) = {[1 0 0 0 -Gload]};
        end
        ForceVec_PCell = cell(num_pat,1);
        for k = 1:num_pat
            p = multiPatSurf_Analysis{k,1};
            q = multiPatSurf_Analysis{k,2};
            mAnalysis = multiPatSurf_Analysis{k,3};
            nAnalysis = multiPatSurf_Analysis{k,4};
            UAnalysis = multiPatSurf_Analysis{k,5};
            VAnalysis = multiPatSurf_Analysis{k,6};
            Ctrlpts = multiPatSurf_Analysis{k,7};
            
            UEle  	= UEle_Cell{k};
            VEle   	= VEle_Cell{k};
            loadType = loadTypeVec(k);
            ploadsMat = ploadsMat_PCell{k};
            ne_u = num_ele_Mat(k,1);
            ne_v = num_ele_Mat(k,2);
            ForceVec = ComputeF(mAnalysis,nAnalysis,p,q,UAnalysis,VAnalysis,UEle,VEle,Ctrlpts,ploadsMat,loadType,ne_u,ne_v);
            ForceVec_PCell(k) = {ForceVec};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Assign displacement degrees of freedom (global indicied of control points in analysis model)
        %%%%%%%%%%%%%%%%%%%%%%% Compute degrees of freedom in analysis model %%%%%%%%%%%%%%%%%%%%%%%%%
        DOF_Analysis = 0;
        Analysis_Patch_DOF_PCell = cell(num_pat,1);
        mAnalysis = multiPatSurf_Analysis{1,3};
        nAnalysis = multiPatSurf_Analysis{1,4};
        DOF_Analysis = DOF_Analysis + (mAnalysis+1) * (nAnalysis+1);
        Analysis_Patch_DOF_PCell(1) = {reshape(1:DOF_Analysis,mAnalysis+1,nAnalysis+1)};
        
        for k = 2:num_pat
            mAnalysis = multiPatSurf_Analysis{k,3};
            nAnalysis = multiPatSurf_Analysis{k,4};
            Analysis_Patch_DOF_PCell(k) = {zeros(mAnalysis+1,nAnalysis+1)};
        end
        
        for k = 1:num_pat - 1
            Patch_DOF = cell2mat(Analysis_Patch_DOF_PCell(k));  % for patch k;
            
            for L = 1:4       % interface repeated dofs
                int = interfaceMat(k,L);
                if (int == 0)
                    continue
                else
                    [row, col] = find(interfaceMat(k+1:end,:) == int);
                    if(isempty(row))
                        continue;
                    else
                        switch L
                            case 1
                                interface_DOF = Patch_DOF(:,1);
                            case 2
                                interface_DOF = Patch_DOF(end,:);
                            case 3
                                interface_DOF = Patch_DOF(:,end);
                            case 4
                                interface_DOF = Patch_DOF(1,:);
                        end
                        
                        row = row + k;
                        Patch_DOF_row = Analysis_Patch_DOF_PCell{row};
                        switch col
                            case 1
                                Patch_DOF_row(:,1) = reshape(interface_DOF,[],1);
                                
                            case 2
                                Patch_DOF_row(end,:) = reshape(interface_DOF,1,[]);
                                
                            case 3
                                Patch_DOF_row(:,end) = reshape(interface_DOF,[],1);
                                
                            case 4
                                Patch_DOF_row(1,:) = reshape(interface_DOF,1,[]);
                                
                        end
                        
                        Analysis_Patch_DOF_PCell(row) = {Patch_DOF_row};
                    end
                end
            end
            
            PreCorPtsNum = (k-1)*4;
            for Cor = 1:4       % corner repeated dofs
                switch Cor
                    case 1
                        DOF_idx = Patch_DOF(1,1);
                    case 2
                        DOF_idx = Patch_DOF(end,1);
                    case 3
                        DOF_idx = Patch_DOF(end,end);
                    case 4
                        DOF_idx = Patch_DOF(1,end);
                end
                
                % find patches that shares the same corner point
                CorPt = CorPtsMat_(PreCorPtsNum+Cor,:);
                diffMat = repmat(CorPt,num_pat*4,1) - CorPtsMat_;
                distMat = sum(diffMat.*diffMat, 2);
                idx_CorPts = find(distMat((PreCorPtsNum+4+1):end) < 1E-3);
                for a = 1:length(idx_CorPts)
                    idx = idx_CorPts(a);
                    idx_pat =ceil(idx/4) + k; % patch k_
                    idx_cor = mod(idx,4);     % point idx in patch k_
                    if(idx_cor == 0)
                        idx_cor = 4;
                    end
                    
                    Patch_DOF_k_ = Analysis_Patch_DOF_PCell{idx_pat};
                    switch idx_cor
                        case 1
                            Patch_DOF_k_(1,1) = DOF_idx;
                        case 2
                            Patch_DOF_k_(end,1) = DOF_idx;
                        case 3
                            Patch_DOF_k_(end,end) = DOF_idx;
                        case 4
                            Patch_DOF_k_(1,end) = DOF_idx;
                    end
                    
                    Analysis_Patch_DOF_PCell(idx_pat) = {Patch_DOF_k_};
                end
                
            end
            
            
            if(k < num_pat)
                Patch_DOF_kPlusOne = Analysis_Patch_DOF_PCell{k+1};  % for patch k;
                interface_idx = find(Patch_DOF_kPlusOne == 0);
                DOF_patch_kPlusOne = (1:numel(interface_idx)) + DOF_Analysis;
                Patch_DOF_kPlusOne(interface_idx) = DOF_patch_kPlusOne;
                Analysis_Patch_DOF_PCell(k+1) = {Patch_DOF_kPlusOne};
                DOF_Analysis = DOF_Analysis + numel(interface_idx);
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Compute global idx for dofs in analysis model -- update ANonZeroIdx_Cell
        %%%%%%%%%%%%%%%%%%%%%%% Compute global idx for dofs in analysis model %%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:num_pat
            ne_u = num_ele_Mat(k,1);
            ne_v = num_ele_Mat(k,2);
            ANonZeroIdx = ANonZeroIdx_PCell{k};
            Analysis_Patch_DOF = Analysis_Patch_DOF_PCell{k};
            for j = 1:ne_v
                for i = 1:ne_u
                    ij = (j-1)*ne_u+i;
                    Idx = ANonZeroIdx(ij,:);
                    Global_Idx = Analysis_Patch_DOF(Idx);
                    ANonZeroIdx(ij,:) = Global_Idx;
                end
            end
            ANonZeroIdx_PCell(k) = {ANonZeroIdx};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Assemble global force vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Assemble the global force vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F_Global = zeros(5*DOF_Analysis,1);
        for k = 1:num_pat
            Patch_DOF_Mat = Analysis_Patch_DOF_PCell{k};
            ForceVec = ForceVec_PCell{k};
            Idx = reshape(((reshape(Patch_DOF_Mat,1,[])-1)*5+(1:5)'),[],1);
            F_Global(Idx) = F_Global(Idx) + ForceVec;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% bdr condition
        switch sType
            case 1
                Patch1_DOF = Analysis_Patch_DOF_PCell{1};
                fixed_DOF = [Patch1_DOF(1,1), Patch1_DOF(1,end), Patch1_DOF(end,end), Patch1_DOF(end,1)];
            case 2
                Patch1_DOF = Analysis_Patch_DOF_PCell{1};
                Patch2_DOF = Analysis_Patch_DOF_PCell{2};
                Patch3_DOF = Analysis_Patch_DOF_PCell{3};
                Patch4_DOF = Analysis_Patch_DOF_PCell{4};
                fixed_DOF = [Patch1_DOF(1,1), Patch2_DOF(end,1), Patch3_DOF(end,end), Patch4_DOF(1,end)];
            case 4
                Patch1_DOF = Analysis_Patch_DOF_PCell{1};
                fixed_DOF = reshape(Patch1_DOF(:,end),1,[]);
                Patch3_DOF = Analysis_Patch_DOF_PCell{3};
                fixed_DOF = [fixed_DOF, reshape(Patch3_DOF(:,1),1,[])];
            case 5
                Patch1_DOF = Analysis_Patch_DOF_PCell{1};
                Patch2_DOF = Analysis_Patch_DOF_PCell{2};
                Patch3_DOF = Analysis_Patch_DOF_PCell{3};
                Patch4_DOF = Analysis_Patch_DOF_PCell{4};
                fixed_DOF = [Patch1_DOF(end,1), Patch1_DOF(1,end), Patch2_DOF(end,1), Patch2_DOF(1,end),Patch3_DOF(end,1), Patch3_DOF(1,end),Patch4_DOF(end,1), Patch4_DOF(1,end)];
            case 6
                fixed_DOF =[];
                for k = 2:4:24
                    Patch_DOF = Analysis_Patch_DOF_PCell{k};
                    fixed_DOF = [fixed_DOF, Patch_DOF(1,:)];
                end
                for k = 3:4:24
                    Patch_DOF = Analysis_Patch_DOF_PCell{k};
                    fixed_DOF = [fixed_DOF, Patch_DOF(1,:)];
                end
            case 3
                Patch1_DOF = cell2mat(Analysis_Patch_DOF_PCell(1));
                Patch2_DOF = cell2mat(Analysis_Patch_DOF_PCell(2));
                Patch3_DOF = cell2mat(Analysis_Patch_DOF_PCell(3));
                Patch4_DOF = cell2mat(Analysis_Patch_DOF_PCell(4));
                fixed_DOF = [Patch1_DOF(1,:), Patch2_DOF(1,:), Patch3_DOF(1,:), Patch4_DOF(1,:)];
        end
        fixed_DOF = unique(fixed_DOF);
        
        %% degree matrix for patches
        DegreeMat = zeros(num_pat, 2);
        for k = 1:num_pat
            p = multiPatSurf_{k,1};
            q = multiPatSurf_{k,2};
            DegreeMat(k,:) = [p q];
        end
        
        %% Design model: Number of design variables in each patch
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Set number of design variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        num_DesignDOF_Mat = zeros(num_pat, 2);
        Mall = max(pat_len_wid_Mat,[],'all');
        for k = 1: num_pat
            len = max(pat_len_wid_Mat(k,1),pat_len_wid_Mat(k,3));
            wid = max(pat_len_wid_Mat(k,2),pat_len_wid_Mat(k,4));
            num_DesignDOF_Mat(k,1) = round(len/Mall*DctrlElenum);
            num_DesignDOF_Mat(k,2) = round(wid/Mall*DctrlElenum);
        end
        
        % uniform design DOFs distribution for patch interfaces according to length
        num_DesignDOF_Mat_ = repmat(num_DesignDOF_Mat,1,2);
        for int = 1:num_interface
            [rows, cols] = find(interfaceMat_ == int);
            if(numel(rows) == 0)
                continue
            else
                ind = sub2ind([num_pat,4],rows,cols);
                max_DesignDOF_num = max(num_DesignDOF_Mat_(ind));
                for ii = 1: numel(rows)
                    row = rows(ii);
                    col = cols(ii);
                    if(mod(col,2) == 1)
                        num_DesignDOF_Mat(row,1) = max_DesignDOF_num;
                    else
                        num_DesignDOF_Mat(row,2) = max_DesignDOF_num;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Refinement for Deisign model:  Weights & UDesign & VDesign
        %%%%%%%%%%%%%%%%%%% Refine the design model through knot insertion %%%%%%%%%%%%%%%%%%%%%
        Design_UV_PCell = cell(num_pat,2);
        Design_Weights_PCell = cell(num_pat,1);
        Design_CtrlPts_Cell = cell(num_pat,1);
        for k = 1:num_pat
            m = cell2mat(multiPatSurf_(k,3));
            n = cell2mat(multiPatSurf_(k,4));
            U = cell2mat(multiPatSurf_(k,5));
            V = cell2mat(multiPatSurf_(k,6));
            Ctrlpts = cell2mat(multiPatSurf_(k,7));
            HCtrlpts = Ctrlpts;
            HCtrlpts(:,:,1:3) = HCtrlpts(:,:,1:3).*HCtrlpts(:,:,4);
            [n1,n2,~]           = size(HCtrlpts);
            % the original num_DesignDOF_Mat(k,:) denote the numbers of elements for design doamin;
            switch sType
                case 4
                    Dctrlnum_u  = num_DesignDOF_Mat(k,1)+3;
                    Dctrlnum_v  = num_DesignDOF_Mat(k,2)+3;
                    tmp = floor((Dctrlnum_u-n1)/2);
                    a = linspace(0,0.5,tmp+2);  a = [a linspace(0.5,1,tmp+2)];
                    tmp = floor((Dctrlnum_v-n2)/2);
                    b = linspace(0,0.5,tmp+2);  b = [b linspace(0.5,1,tmp+2)];
                otherwise
                    ne_s_more  = num_DesignDOF_Mat(k,1) - (numel(unique(U))-1);
                    ne_t_more  = num_DesignDOF_Mat(k,2) - (numel(unique(V))-1);
                    a = setdiff(linspace(U(1),U(end),ne_s_more+2),U);
                    b = setdiff(linspace(V(1),V(end),ne_t_more+2),V);
            end
            
            a = setdiff(a,U); b = setdiff(b,V);
            %             a = [];b = [];
            UDesign    	= sort([U,a]);
            VDesign   	= sort([V,b]);
            Dctrlnum_u  = (m+1) + length(UDesign) - length(U);
            Dctrlnum_v  = (n+1) + length(VDesign) - length(V);
            num_DesignDOF_Mat(k,:) = [Dctrlnum_u, Dctrlnum_v];
            % the num_DesignDOF_Mat(k,:) are used to represent the number of deisgn variables for design model here.(dofs)
            
            [~, ~, DHCtrlpts]    = KnotInsert(a,b,U,V,HCtrlpts);
            DCtrlpts = DHCtrlpts(:,:,1:3)./DHCtrlpts(:,:,4);
            Weights = reshape(DHCtrlpts(:,:,4), Dctrlnum_u,Dctrlnum_v);
            Design_CtrlPts_Cell(k) = {DCtrlpts};
            
            Design_UV_PCell(k,:) = {UDesign,VDesign};
            Design_Weights_PCell(k) = {Weights};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% DOF in Design model
        %%%%%%%%%%%%%%%%%%%%%%% Set degree of freedom in design model %%%%%%%%%%%%%%%%%%%%%%%%%%
        DOF_Design = 0;
        Design_Patch_DOF_PCell = cell(num_pat,1);
        mDesign = num_DesignDOF_Mat(1,1)-1;
        nDesign = num_DesignDOF_Mat(1,2)-1;
        DOF_Design = DOF_Design + (mDesign+1) * (nDesign+1);
        Design_Patch_DOF_PCell(1) = {reshape(1:DOF_Design,mDesign+1,nDesign+1)};
        
        for k = 2:num_pat
            Design_Patch_DOF_PCell(k) = {zeros(num_DesignDOF_Mat(k,1),num_DesignDOF_Mat(k,2))};
        end
        
        for k = 1:num_pat - 1
            Design_Patch_DOF = cell2mat(Design_Patch_DOF_PCell(k));  % for patch k;
            
            for L = 1:4     % interface repeated dofs
                int = interfaceMat(k,L);
                if (int == 0)
                    continue
                else
                    [row, col] = find(interfaceMat(k+1:end,:) == int);
                    if(isempty(row))
                        continue;
                    else
                        
                        switch L
                            case 1
                                interface_DOF = Design_Patch_DOF(1:end,1);
                            case 2
                                interface_DOF = Design_Patch_DOF(end,1:end);
                            case 3
                                interface_DOF = Design_Patch_DOF(1:end,end);
                            case 4
                                interface_DOF = Design_Patch_DOF(1,1:end);
                        end
                        
                        row = row + k;
                        Patch_DOF_row = cell2mat(Design_Patch_DOF_PCell(row));
                        switch col
                            case 1
                                Patch_DOF_row(1:end,1) = reshape(interface_DOF,[],1);
                            case 2
                                Patch_DOF_row(end,1:end) = reshape(interface_DOF,1,[]);
                            case 3
                                Patch_DOF_row(1:end,end) = reshape(interface_DOF,[],1);
                            case 4
                                Patch_DOF_row(1,1:end) = reshape(interface_DOF,1,[]);
                        end
                        
                        Design_Patch_DOF_PCell(row) = {Patch_DOF_row};
                    end
                end
            end
            
            PreCorPtsNum = (k-1)*4;
            for Cor = 1:4   % corner repeated dofs
                switch Cor
                    case 1
                        DOF_idx = Design_Patch_DOF(1,1);
                    case 2
                        DOF_idx = Design_Patch_DOF(end,1);
                    case 3
                        DOF_idx = Design_Patch_DOF(end,end);
                    case 4
                        DOF_idx = Design_Patch_DOF(1,end);
                end
                
                % find patches that shares the same corner point
                CorPt = CorPtsMat_(PreCorPtsNum+Cor,:);
                diffMat = repmat(CorPt,num_pat*4,1) - CorPtsMat_;
                distMat = sum(diffMat.*diffMat, 2);
                idx_CorPts = find(distMat((PreCorPtsNum+4+1):end) < 1E-3);
                for a = 1:length(idx_CorPts)
                    idx = idx_CorPts(a);
                    idx_pat =ceil(idx/4) + k; % patch k_
                    idx_cor = mod(idx,4);     % point idx in patch k_
                    if(idx_cor == 0)
                        idx_cor = 4;
                    end
                    
                    Patch_DOF_k_ = Design_Patch_DOF_PCell{idx_pat};
                    switch idx_cor
                        case 1
                            Patch_DOF_k_(1,1) = DOF_idx;
                        case 2
                            Patch_DOF_k_(end,1) = DOF_idx;
                        case 3
                            Patch_DOF_k_(end,end) = DOF_idx;
                        case 4
                            Patch_DOF_k_(1,end) = DOF_idx;
                    end
                    
                    Design_Patch_DOF_PCell(idx_pat) = {Patch_DOF_k_};
                end
                
            end
            
            if(k < num_pat)
                Design_Patch_DOF_kPlusOne = cell2mat(Design_Patch_DOF_PCell(k+1));  % for patch k;
                interface_idx = find(Design_Patch_DOF_kPlusOne == 0);
                Design_DOF_patch_kPlusOne = (1:numel(interface_idx)) + DOF_Design;
                Design_Patch_DOF_kPlusOne(interface_idx) = Design_DOF_patch_kPlusOne;
                Design_Patch_DOF_PCell(k+1) = {Design_Patch_DOF_kPlusOne};
                DOF_Design = DOF_Design + numel(interface_idx);
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Filter Type: LVC / NULL
        EleInfoMat = [];                            % Record Element Information: (ele_pat, ele_prow, ele_pcol);
        Neighbors_ECell = [];                       % Record the indicies of neighbor-elements for each element.
        EleIdx_PCell = cell(num_pat,1);             % Record the global ele-idx for elements of each patch
        
        switch myFilter
            case 'NULL'
                EleIdx_PCell = {};
            case 'LocalVolumeConstraint'
                % Information for Local Volume Constraint
                [EleInfoMat, EleIdx_PCell, Neighbors_ECell] = Find_EleNeighbors_LVC(multiPatSurf,num_pat,interfaceMat,num_ele_Mat,pat_len_wid_Mat,UEle_Cell,VEle_Cell);
        end
        
        %% basis in Gauss points
        NonZeroBasis_GuassPt_PCell = cell(num_pat,1);
        DNonZeroIdx_GuassPt_PCell = cell(num_pat,1);
        for k = 1:num_pat
            p = cell2mat(multiPatSurf_(k,1));
            q = cell2mat(multiPatSurf_(k,2));
            Dctrlnum_u = num_DesignDOF_Mat(k,1);
            Dctrlnum_v = num_DesignDOF_Mat(k,2);
            UDesign = cell2mat(Design_UV_PCell(k,1));
            VDesign = cell2mat(Design_UV_PCell(k,2));
            UEle = cell2mat(UEle_Cell(k));
            VEle = cell2mat(VEle_Cell(k));
            Weights = cell2mat(Design_Weights_PCell(k));
            ne_u = num_ele_Mat(k,1);
            ne_v = num_ele_Mat(k,2);
            [NonZeroBasis_GuassPt_ECell, DNonZeroIdx_ECell] = ComputeNonZeroBasisAndIdx(p,q,Dctrlnum_u,Dctrlnum_v,UDesign,VDesign,UEle,VEle,Weights,ne_u,ne_v);
            NonZeroBasis_GuassPt_PCell(k) = {NonZeroBasis_GuassPt_ECell};
            
            Design_Patch_DOF = Design_Patch_DOF_PCell{k};
            for j = 1:ne_v
                for i = 1:ne_u
                    ij = (j-1)*ne_u+i;
                    DNonZeroIdx = DNonZeroIdx_ECell{ij};
                    for kgp = 1:ng_u*ng_v
                        Idx = DNonZeroIdx(kgp,:);
                        Global_Idx = Design_Patch_DOF(Idx);
                        DNonZeroIdx_ECell{ij}(kgp,:) = Global_Idx;
                    end
                end
            end
            DNonZeroIdx_GuassPt_PCell(k) = {DNonZeroIdx_ECell};
        end
        
        %% topology optimization using mma
        VOL = 0;
        for k = 1:num_pat
            Ve_GuassPt_All = StiffnessVolWeight_PCell{k,2};
            Weight_GuassPt_All = StiffnessVolWeight_PCell{k,3};
            VOL = VOL + sum(Ve_GuassPt_All.*Weight_GuassPt_All,'all');
        end
        
        switch myFilter
            case 'NULL'
                volfrac = 0.4;
                VolumeCons = volfrac * VOL;
            case 'LocalVolumeConstraint'
                volfrac = 0.6;
                VolumeCons = volfrac;
        end
        
        nxval = DOF_Design;
        xval = ones(nxval,1);
        
        num_constraints = 1;
        fval = zeros(num_constraints,1);
        
        a0 = 1;
        a = zeros(num_constraints,1);
        c = ones(num_constraints,1) * 1000;
        d = zeros(num_constraints,1);
        
        xold1 = xval;
        xold2 = xval;
        xmin  = 0.001 * ones(nxval,1);
        xmax  = ones(nxval,1);
        low   = xmin;
        upp   = xmax;
        
        %% iteration for optimization process
        beta = 2;
        loop = 0;
        betascale = 0;
        Data = zeros(200,4);
        
        clear NonZeroBasis_GuassPt_ECell DNonZeroIdx_ECell ANonZeroIdx Analysis_Patch_DOF Design_CtrlPts_Cell Design_Patch_DOF Design_Patch_DOF_kPlusOne ForceVec ForceVec_Cell
        clear DOF_patch_kPlusOne Patch_DOF Patch_DOF_kPlusOne Patch_DOF_Mat PatchK ploadsMat_PCell Ve_GuassPt_All  Ke_GuassPt_All Weight_GuassPt_All Analysis_Patch_DOF_PCell
        F_Global = sparse(F_Global);
        while (loop < 200)
            tstart = tic;
            loop = loop + 1;
            betascale = betascale + 1;
            if(beta < 64 && betascale == 25)
                beta = beta * 2;
                betascale = 0;
            end
            DensCtrlpts = xval;
            [compl,DersCompl_Vec,vol,DersVol_Vec,GlobalVol] = IGAInfo(DOF_Analysis,DOF_Design,num_pat,DegreeMat,num_ele_Mat,DensCtrlpts,NonZeroBasis_GuassPt_PCell,DNonZeroIdx_GuassPt_PCell,...
                ANonZeroIdx_PCell,StiffnessVolWeight_PCell,F_Global,fixed_DOF,myFilter,EleInfoMat,EleIdx_PCell,Neighbors_ECell);
            
            %% test accuracy of Derivatives by FDM
            %         % 1.DerVol
            %         [compl,DersCompl_Vec,vol,DersVol_Vec] = IGAInfo(DOF_Analysis,DOF_Design,num_pat,DegreeMat,num_ele_Mat,DensCtrlpts,NonZeroBasis_GuassPt_PCell,DNonZeroIdx_GuassPt_PCell,...
            %                                                         ANonZeroIdx_PCell,StiffnessVolWeight_PCell,F_Global,fixed_DOF,myFilter,EleInfoMat,EleIdx_PCell,Neighbors_ECell);
            %         partial_rho_C = DersCompl_Vec(3);
            %         partial_rho_V = DersVol_Vec(3);
            %
            %         % 2.DerVol FDM
            %         mn = numel(DensCtrlpts);
            %         FDM_C = zeros(mn,1);
            %         FDM_V = zeros(mn,1);
            %         for i = 1:numel(DensCtrlpts)
            %             DensCtrlpts_ = DensCtrlpts;
            %             DensCtrlpts_(i) = DensCtrlpts(i) + (1e-4);
            %             [compl_,DersCompl_Vec,vol_,DersVol_Vec] = IGAInfo(DOF_Analysis,DOF_Design,num_pat,DegreeMat,num_ele_Mat,DensCtrlpts_,NonZeroBasis_GuassPt_PCell,DNonZeroIdx_GuassPt_PCell,...
            %                                                               ANonZeroIdx_PCell,StiffnessVolWeight_PCell,F_Global,fixed_DOF,myFilter,EleInfoMat,EleIdx_PCell,Neighbors_ECell);
            %             partial_rho_C_ = (compl_ - compl)/(1e-4);
            %             partial_rho_V_ = (vol_ - vol)/(1e-4);
            %             FDM_C(i) = partial_rho_C_;
            %             FDM_V(i) = partial_rho_V_;
            %         end
            %         AAA = DersCompl_Vec - FDM_C;
            %         BBB = DersVol_Vec - FDM_V;
            
            %% solver
            cScale = 100 / compl;
            f0val = 100;
            df0dx = DersCompl_Vec * cScale;
            fval  = vol/VolumeCons - 1;
            dfdx  = DersVol_Vec' / VolumeCons;
            
            [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(num_constraints,nxval,loop,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
            xold2 = xold1;
            xold1 = xval;
            xval = 0.5*xold1 + 0.5*xmma;
            
            time = toc(tstart);
            fprintf('Iter: %d, Obj: %f, Con: %f, Time: %f\n',loop,compl,fval,time);
            Data(loop,:) = [loop,compl,fval,time];
        end
        
        %% Draw density of mid-surface
        figure(kij)
        set(figure(kij),'Position',[300 300 350 280]);
        multiple = ceil(50/min(num_ele_Mat,[],'all'));
        for k = 1:num_pat
            Design_Patch_DOF = cell2mat(Design_Patch_DOF_PCell(k));
            DensCtrlpts = xval(Design_Patch_DOF);
            [Dctrlnum_s,Dctrlnum_t] = size(DensCtrlpts);
            UDesign = cell2mat(Design_UV_PCell(k,1));
            VDesign = cell2mat(Design_UV_PCell(k,2));
            
            p = cell2mat(multiPatSurf_Analysis(k,1));
            q = cell2mat(multiPatSurf_Analysis(k,2));
            m = cell2mat(multiPatSurf_Analysis(k,3));
            n = cell2mat(multiPatSurf_Analysis(k,4));
            U_Analysis = cell2mat(multiPatSurf_Analysis(k,5));
            V_Analysis = cell2mat(multiPatSurf_Analysis(k,6));
            Ctrlpts = cell2mat(multiPatSurf_Analysis(k,7));
            
            
            Nx = multiple*num_ele_Mat(k,1)+1;
            Ny = multiple*num_ele_Mat(k,2)+1;
            u = linspace(U(1),U(end),Nx);
            v = linspace(V(1),V(end),Ny);
            pnum = Nx*Ny;
            surf = zeros(pnum,3);
            Density = zeros(Nx,Ny);
            
            for j = 1:Ny
                for i = 1:Nx
                    ii = (j-1)*Nx+i;
                    [Val,IdxU,IdxV] = NonZeroBasis(m, n, p, q, u(i), v(j), U_Analysis, V_Analysis, reshape(Ctrlpts(:,:,4),m+1,n+1));
                    P               = sum(Ctrlpts(IdxU+1,IdxV+1,1:3).*Val,[1,2]);
                    surf(ii,:)       = reshape(P,1,3);
                    
                    [Val,IdxU,IdxV] = NonZeroBasis(Dctrlnum_s-1, Dctrlnum_t-1,p, q, u(i), v(j), UDesign, VDesign, cell2mat(Design_Weights_PCell(k)));
                    den             = sum(DensCtrlpts(IdxU+1,IdxV+1).*Val,'all');
                    Density(i,j)    = Heaviside(den,beta,0.5);
                end
            end
            
            figure(kij)
            Px = reshape(surf(:,1),Nx,Ny);
            Py = reshape(surf(:,2),Nx,Ny);
            Pz = reshape(surf(:,3),Nx,Ny);
            Cden = zeros((Nx-1)*(Ny-1),4);
            lx = zeros((Nx-1)*(Ny-1),4);
            ly = zeros((Nx-1)*(Ny-1),4);
            lz = zeros((Nx-1)*(Ny-1),4);
            for i = 1:Nx-1
                for j = 1:Ny-1
                    lx((i-1)*(Ny-1)+j,1:4)=[Px(i,j) Px(i+1,j) Px(i+1,j+1) Px(i,j+1)];
                    ly((i-1)*(Ny-1)+j,1:4)=[Py(i,j) Py(i+1,j) Py(i+1,j+1) Py(i,j+1)];
                    lz((i-1)*(Ny-1)+j,1:4)=[Pz(i,j) Pz(i+1,j) Pz(i+1,j+1) Pz(i,j+1)];
                    Cden((i-1)*(Ny-1)+j,1:4)=[Density(i,j) Density(i+1,j) Density(i+1,j+1) Density(i,j+1)];
                end
            end
            
            for i=1:size(lx,1)
                hold on
                fill3(lx(i,:)',ly(i,:)',lz(i,:)', Cden(i,:),'LineStyle','none');
            end
            
        end
        
        figure(kij)
        set(gca,'CLim',[0,1]);
        colorbar()
        txt = ['Obj=' num2str(compl), '  Con=' num2str(fval)];
        xlabel(txt)
        axis equal
        title(sprintf("Cpts = %d %d", num_DesignDOF_Mat(1,1), num_DesignDOF_Mat(1,2)));
        text(mean(xlim()), min(ylim()), txt, 'HorizontalAlignment', 'center');
        g = gray;
        g = flipud(g);
        colormap(g);
        colorbar('off') ;
        axis off;
        
        %% save and clear
        currentDir = pwd;
        subfolder = ['s',num2str(sType)];    % subfolder
        if ~exist(subfolder, 'dir')             % create the subfolder if does not exist
            mkdir(subfolder);
        end
        
        if(num_pat < 3)
            path = fullfile(currentDir, subfolder,['DensCtrlpts = ', num2str(num_DesignDOF_Mat(1,1)),' ', num2str(num_DesignDOF_Mat(1,2)) ,' num_Ele= ', num2str(num_ele)]);
        else
            path = fullfile(currentDir, subfolder,['DensCtrlpts = ', num2str(num_DesignDOF_Mat(1,1)),' ', num2str(num_DesignDOF_Mat(3,1)) ,' ', num2str(num_DesignDOF_Mat(1,2)) ,' ', num2str(num_DesignDOF_Mat(2,2)),' num_Ele= ', num2str(num_ele)]);
        end
        
        filename = [path,'_DensCpts'];
        save(path, 'xval', 'Design_UV_PCell','Design_Weights_PCell','Design_Patch_DOF_PCell','num_ele_Mat','interfaceMat','Data','multiPatSurf_');
        saveas(figure(kij),filename,'jpg');
        saveas(figure(kij),filename,'fig');
        close(figure(kij));
        
        filename = [path,'_Data'];
        xlswrite(sprintf("%s.xls",filename), Data);
        
        
        filename = [path,'_NumDesignDOF'];
        xlswrite(sprintf("%s .xls",filename), num_DesignDOF_Mat);
        
        filename = [path,'_Compliance'];
        figure(15+kij)
        plot(Data(:,1), Data(:,2),'b-');
        title("compliance-iteration");
        xlabel('iteration');
        ylabel('compliance');
        saveas(figure(15+kij),filename,'jpg');
        close(figure(15+kij));
        
        filename = [path,'_Volume'];
        figure(30+kij)
        plot(Data(:,1), Data(:,3),'b-');
        title("volume-iteration");
        xlabel('iteration');
        ylabel('volume');
        saveas(figure(30+kij),filename,'jpg');
        close(figure(30+kij));
        
        temp = who; % get all variable names in the current workspace
        vars_to_keep = {'kij','myFilter'};
        vars_to_clear = setdiff(temp, vars_to_keep);
        clear(vars_to_clear{:});
        clc;
        
    
end

end

function ANonZeroIdx = ComputeNonZeroBasis(p,q,mAnalysis,nAnalysis,UAnalysis,VAnalysis,UEle,VEle,ne_u,ne_v)
ANonZeroIdx = zeros(ne_u*ne_v,(p+1)*(q+1));
for j = 1:ne_v
    t = (VEle(j)+VEle(j+1))/2;
    for i = 1:ne_u
        s = (UEle(i)+UEle(i+1))/2;
        ij = (j-1)*ne_u+i;
        IdxU = NonZeroBasisIdx(mAnalysis, p, s, UAnalysis);
        IdxV = NonZeroBasisIdx(nAnalysis, q, t, VAnalysis);
        [J, I] = meshgrid(IdxV, IdxU);
        ANonZeroIdx(ij,:) = sub2ind([mAnalysis+1,nAnalysis+1], I(:)+1, J(:)+1);
    end
end
end

function [NonZeroBasis_GuassPt_ECell, DNonZeroIdx_ECell] = ComputeNonZeroBasisAndIdx(p,q,Dctrlnum_u,Dctrlnum_v,UDesign,VDesign,UEle,VEle,Weights,ne_u,ne_v)
global ng_u ng_v;  % Emid --> E_GuassPts
NonZeroBasis_GuassPt_ECell = cell(ne_u*ne_v,1);
DNonZeroIdx_ECell      = cell(ne_u*ne_v,1);

NonZeroBasis_GuassPt = zeros(ng_u*ng_v,(p+1)*(q+1));
DNonZeroIdx      = zeros(ng_u*ng_v,(p+1)*(q+1));

for j = 1:ne_v
    [s2, ~] = GaussInt(VEle(j),VEle(j+1),ng_v);
    for i = 1:ne_u
        [s1, ~] = GaussInt(UEle(i),UEle(i+1),ng_u);
        ij = (j-1)*ne_u+i;
        
        for kgp = 1:ng_u*ng_v
            [tmp1, tmp2] = meshgrid(s1,s2);
            uv = [reshape(tmp1,[],1),reshape(tmp2,[],1)];
            [Val, IdxU, IdxV] = NonZeroBasis(Dctrlnum_u-1, Dctrlnum_v-1, p, q, uv(kgp,1), uv(kgp,2), UDesign, VDesign, Weights);
            [J, I] = meshgrid(IdxV, IdxU);
            Idx = sub2ind([Dctrlnum_u, Dctrlnum_v], I(:)+1, J(:)+1);
            NonZeroBasis_GuassPt(kgp,:) = reshape(Val,[],1);
            DNonZeroIdx(kgp,:) = Idx;
        end
        
        NonZeroBasis_GuassPt_ECell(ij) = {NonZeroBasis_GuassPt};
        DNonZeroIdx_ECell(ij) = {DNonZeroIdx};
    end
end
end

function [EleInfoMat, EleIdx_PCell, Neighbors_ECell] = Find_EleNeighbors_LVC(multiPatSurf,num_pat,interfaceMat,num_ele_Mat,pat_len_wid_Mat,UEle_Cell,VEle_Cell)
% element neighbors for local volume constraint
global filter_radius;
Radius = filter_radius * sum(pat_len_wid_Mat,'all')/sum(num_ele_Mat(:,1)+num_ele_Mat(:,2))/2;
EleIdx_PCell = cell(num_pat,1); % Record the global ele-idx for elements of each patch
ElePos_PCell = cell(num_pat,1); % Record the center position for elementsof each patch
sum_num_ele = 0;
for k = 1:num_pat
    ne_u = num_ele_Mat(k,1);
    ne_v = num_ele_Mat(k,2);
    EleIdx_PCell(k) = {sum_num_ele+reshape(1:(ne_u*ne_v),ne_u,ne_v)};
    sum_num_ele = sum_num_ele + ne_u*ne_v;
    
    ElePos = zeros(ne_u*ne_v,3);
    UEle = UEle_Cell{k};
    VEle = VEle_Cell{k};
    p = multiPatSurf{k,1};
    q = multiPatSurf{k,2};
    m = multiPatSurf{k,3};
    n = multiPatSurf{k,4};
    U = multiPatSurf{k,5};
    V = multiPatSurf{k,6};
    Ctrlpts = multiPatSurf{k,7};
    Weights = reshape(Ctrlpts(:,:,4),m+1,n+1);
    for j = 1:ne_v
        v = (VEle(j)+VEle(j+1))/2;
        for i = 1:ne_u
            u = (UEle(i)+UEle(i+1)) / 2;
            [IdxU, IdxV, N, ~, ~, ~, ~, ~] = NrbSurDer2(m,n,p,q,u,v,U,V,Weights);
            PIJ = Ctrlpts(IdxU+1,IdxV+1,1:3);
            ElePos((j-1)*ne_u+i,:) = reshape(sum(PIJ.*N, [1,2]),1,3);
        end
    end
    
    ElePos_PCell(k) = {ElePos};
end

% sum_ele = sum_num_ele;
sum_ele = sum(num_ele_Mat(:,1).*num_ele_Mat(:,2)); % total number of elements
EleInfoMat = zeros(sum_ele, 3); % Record Element Information: (ele_pat, ele_prow, ele_pcol);
Neighbors_ECell = cell(sum_ele, 1); % Record the indicies of neighbor-elements for each element.
ele_idx = 0;
for k = 1:num_pat
    ne_u = num_ele_Mat(k,1);
    ne_v = num_ele_Mat(k,2);
    EleIdx_Mat = EleIdx_PCell{k};
    ElePos_Mat = ElePos_PCell{k};
    for j = 1:ne_v
        for i = 1:ne_u
            ij = (j-1)*ne_u+i;
            Px = ElePos_Mat(ij,1);
            Py = ElePos_Mat(ij,2);
            Pz = ElePos_Mat(ij,3);
            ele_idx = ele_idx + 1;
            ele_neighbors = zeros(sum_ele,1);
            EleInfoMat(ele_idx, :) = [k,i,j];
            
            % find neighbor elements within the same patch
            dist = ( (ElePos_Mat(:,1)-Px).^2 + (ElePos_Mat(:,2)-Py).^2 + (ElePos_Mat(:,3)-Pz).^2 );
            ele_indicies = find(dist <= Radius^2);
            num_neighbors = numel(ele_indicies);
            ele_neighbors(1:num_neighbors) = EleIdx_Mat(ele_indicies);
            
            
            % find neighbor elements on adjacent patchs through interface
            for edge = 1:4
                int = interfaceMat(k,edge);
                if(int == 0)
                    continue;
                else
                    [rows, cols] = find(interfaceMat == int);
                    for iii = 1:numel(rows) % 1,2
                        pat_kk = rows(iii); % interface patch
                        edge_kk = cols(iii);
                        if(pat_kk == k && edge_kk == edge)
                            continue;
                        else
                            ElePos_Mat_kk = ElePos_PCell{pat_kk};
                            EleIdx_Mat_kk = EleIdx_PCell{pat_kk};
                            dist = ( (ElePos_Mat_kk(:,1)-Px).^2 + (ElePos_Mat_kk(:,2)-Py).^2 + (ElePos_Mat_kk(:,3)-Pz).^2 );
                            ele_indicies = find(dist <= Radius^2);
                            num_neighbors_kk = numel(ele_indicies);
                            ele_neighbors(num_neighbors+(1:num_neighbors_kk)) = EleIdx_Mat_kk(ele_indicies);
                            num_neighbors = num_neighbors + num_neighbors_kk;
                        end
                    end
                end
                
            end
            
            Neighbors_ECell(ele_idx) = {ele_neighbors(1:num_neighbors)};
            
        end
    end
end
end

