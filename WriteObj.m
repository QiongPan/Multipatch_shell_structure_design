function WriteObj()

% read surface and density information
filename = '.\s7\gamma=32\vol = 0.6\LVC_Filter = 6\DensCtrlpts = 33 33 30 30 num_Ele= 80';
sType = 7;

filename_obj =  [filename,'.obj'];
filename_stl = [filename,'.stl'];
filename_mat = [filename,'.mat'];
load(filename_mat); 
% xval = table2array(xval); 

%% multipatch surface
beta = 64;
multiPatSurf = NurbsSurface(sType); 
if sType == 6  % the degree for original model and design model is different
    multiPatSurf = multiPatSurf_;
end

num_pat     = size(multiPatSurf,1);
if num_pat >= 8 
    multiple = ceil(100/min(num_ele_Mat,[],'all'));
elseif num_pat >= 4
    multiple = ceil(200/min(num_ele_Mat,[],'all'));
else
    multiple = ceil(400/min(num_ele_Mat,[],'all'));
end
num_ele_Mat = multiple*num_ele_Mat;
    
% corner point np*4*4
CorPtsMat = zeros(num_pat, 4, 3); 
for k = 1:num_pat
    CPs = multiPatSurf{k,7};
    m_k = multiPatSurf{k,3};
    n_k = multiPatSurf{k,4};
    CorPtsMat(k,:,:) = reshape([CPs(1,1,1:3); CPs(m_k+1,1,1:3); CPs(m_k+1,n_k+1,1:3); CPs(1,n_k+1,1:3)],4,3);
end
CorPtsMat_ = zeros(num_pat*4,3);
CorPtsMat_(:,1) = reshape(CorPtsMat(:,:,1)',num_pat*4,1);
CorPtsMat_(:,2) = reshape(CorPtsMat(:,:,2)',num_pat*4,1);
CorPtsMat_(:,3) = reshape(CorPtsMat(:,:,3)',num_pat*4,1);


%% mesh -- set idx for verticies in each patch
num_vertex = 0;
num_tri = 0;
Nx = num_ele_Mat(1,1)+1;
Ny = num_ele_Mat(1,2)+1;
GridVertIdx_PCell = cell(num_pat,1);
CenterVertIdx_PCell = cell(num_pat,1);
num_vertex = num_vertex + Nx*Ny;
GridVertIdx_PCell(1) = {reshape(1:num_vertex,Nx,Ny)};
CenterVertIdx_PCell(1) = {reshape(1:(Nx-1)*(Ny-1),Nx-1,Ny-1)+num_vertex};
num_vertex = num_vertex +(Nx-1)*(Ny-1);
num_tri = num_tri + 4*(Nx-1)*(Ny-1);
for k = 2:num_pat
    Nx = num_ele_Mat(k,1)+1;
    Ny = num_ele_Mat(k,2)+1;
    GridVertIdx_PCell(k) = {zeros(Nx,Ny)};
    CenterVertIdx_PCell(k) = {reshape(1:(Nx-1)*(Ny-1),Nx-1,Ny-1)};
end

for k = 1:num_pat - 1
    Patch_GridVertIdx = GridVertIdx_PCell{k};
    
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
                        interface_GVIdx = Patch_GridVertIdx(:,1);
                    case 2
                        interface_GVIdx = Patch_GridVertIdx(end,:);
                    case 3
                        interface_GVIdx = Patch_GridVertIdx(:,end);
                    case 4
                        interface_GVIdx = Patch_GridVertIdx(1,:);
                end
            
                row = row + k;
                interface_GVIdx_row = GridVertIdx_PCell{row};
                switch col
                    case 1
                        interface_GVIdx_row(:,1) = reshape(interface_GVIdx,[],1);
                    case 2
                        interface_GVIdx_row(end,:) = reshape(interface_GVIdx,1,[]);
                    case 3
                        interface_GVIdx_row(:,end) = reshape(interface_GVIdx,[],1);
                    case 4
                        interface_GVIdx_row(1,:) = reshape(interface_GVIdx,1,[]);
                end
                GridVertIdx_PCell(row) = {interface_GVIdx_row};
            end
        end
    end
    
    % find patches that shares the same corner point
    PreCorPtsNum = (k-1)*4;
    for Cor = 1:4       % corner repeated dofs
        switch Cor
            case 1
                Vertex_idx = Patch_GridVertIdx(1,1);
            case 2
                Vertex_idx = Patch_GridVertIdx(end,1);
            case 3
                Vertex_idx = Patch_GridVertIdx(end,end);
            case 4
                Vertex_idx = Patch_GridVertIdx(1,end);
        end
                
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
                    
            Patch_GridVertIdx_k_ = GridVertIdx_PCell{idx_pat};
            switch idx_cor
                case 1
                    Patch_GridVertIdx_k_(1,1) = Vertex_idx;
                case 2
                    Patch_GridVertIdx_k_(end,1) = Vertex_idx;
                case 3
                    Patch_GridVertIdx_k_(end,end) = Vertex_idx;
                case 4
                    Patch_GridVertIdx_k_(1,end) = Vertex_idx;
            end
                    
            GridVertIdx_PCell(idx_pat) = {Patch_GridVertIdx_k_};
        end
                
    end
            
    if(k < num_pat)
        Patch_GridVertIdx_kPlusOne = GridVertIdx_PCell{k+1};  % for patch k+1;
        interface_idx = find(Patch_GridVertIdx_kPlusOne == 0);
        Patch_GridVertIdx_kPlusOne(interface_idx) = (1:numel(interface_idx)) + num_vertex;
        GridVertIdx_PCell(k+1) = {Patch_GridVertIdx_kPlusOne};
        num_vertex = num_vertex + numel(interface_idx);
        
        Nx = num_ele_Mat(k+1,1)+1;
        Ny = num_ele_Mat(k+1,2)+1;
        CenterVertIdx_PCell(k+1) = {CenterVertIdx_PCell{k+1}+num_vertex};
        num_vertex = num_vertex + (Nx-1)*(Ny-1);
        num_tri = num_tri + 4*(Nx-1)*(Ny-1);
    end
            
end


%% mesh -- Verticies coordinate & triangular_centerpoint_density & faces
Vertices = zeros(num_vertex,3);
density = zeros(num_tri,1);
Faces = zeros(num_tri,3);
tri_idx = 0;
face_idx = 0;
for k = 1:num_pat
    Design_Patch_DOF = cell2mat(Design_Patch_DOF_PCell(k));
    DensCtrlpts = xval(Design_Patch_DOF);
    [Dctrlnum_s,Dctrlnum_t] = size(DensCtrlpts);
    UDesign = cell2mat(Design_UV_Cell(k,1));
    VDesign = cell2mat(Design_UV_Cell(k,2));
    
    Patch_GridVertIdx = GridVertIdx_PCell{k};
    Patch_CenterVertIdx = CenterVertIdx_PCell{k};
            
    p = multiPatSurf{k,1};
    q = multiPatSurf{k,2};
    m = multiPatSurf{k,3};
    n = multiPatSurf{k,4};
    U = multiPatSurf{k,5};
    V = multiPatSurf{k,6};
    Ctrlpts = multiPatSurf{k,7};
            
    Nx = num_ele_Mat(k,1)+1;
    Ny = num_ele_Mat(k,2)+1;
    uvec = linspace(U(1),U(end),Nx);
    vvec = linspace(V(1),V(end),Ny);
    
    %vertex
    for j = 1:Ny
        for i = 1:Nx
            vertex_idx = Patch_GridVertIdx(i,j);
            [Val,IdxU,IdxV] = NonZeroBasis(m, n, p, q, uvec(i), vvec(j), U, V, reshape(Ctrlpts(:,:,4),m+1,n+1));
            P              = sum(Ctrlpts(IdxU+1,IdxV+1,1:3).*Val,[1,2]);
            Vertices(vertex_idx,:) = reshape(P,1,3);
        end
    end
    
    % vertex and density
    for j = 1:Ny-1
        for i = 1:Nx-1
            vertex_idx = Patch_CenterVertIdx(i,j);
            [Val,IdxU,IdxV] = NonZeroBasis(m, n, p, q, 0.5*uvec(i)+0.5*uvec(i+1), 0.5*vvec(j)+0.5*vvec(j+1), U, V, reshape(Ctrlpts(:,:,4),m+1,n+1));
            P              = sum(Ctrlpts(IdxU+1,IdxV+1,1:3).*Val,[1,2]);
            Vertices(vertex_idx,:) = reshape(P,1,3);
            
            [Val,IdxU,IdxV] = NonZeroBasis(Dctrlnum_s-1, Dctrlnum_t-1, p, q, 0.5*uvec(i)+0.5*uvec(i+1), 0.75*vvec(j)+0.25*vvec(j+1), UDesign, VDesign, Design_Weights_Cell{k});
            den_1              = sum(DensCtrlpts(IdxU+1,IdxV+1).*Val,'all');
            [Val,IdxU,IdxV] = NonZeroBasis(Dctrlnum_s-1, Dctrlnum_t-1, p, q, 0.25*uvec(i)+0.75*uvec(i+1), 0.5*vvec(j)+0.5*vvec(j+1), UDesign, VDesign, Design_Weights_Cell{k});
            den_2               = sum(DensCtrlpts(IdxU+1,IdxV+1).*Val,'all');
            [Val,IdxU,IdxV] = NonZeroBasis(Dctrlnum_s-1, Dctrlnum_t-1, p, q, 0.5*uvec(i)+0.5*uvec(i+1), 0.25*vvec(j)+0.75*vvec(j+1), UDesign, VDesign, Design_Weights_Cell{k});
            den_3               = sum(DensCtrlpts(IdxU+1,IdxV+1).*Val,'all');
            [Val,IdxU,IdxV] = NonZeroBasis(Dctrlnum_s-1, Dctrlnum_t-1, p, q, 0.75*uvec(i)+0.25*uvec(i+1), 0.5*vvec(j)+0.5*vvec(j+1), UDesign, VDesign, Design_Weights_Cell{k});
            den_4               = sum(DensCtrlpts(IdxU+1,IdxV+1).*Val,'all');
            density(tri_idx + (1:4)) = [Heaviside(den_1,beta,0.5); Heaviside(den_2,beta,0.5); Heaviside(den_3,beta,0.5); Heaviside(den_4,beta,0.5)];
            tri_idx = tri_idx + 4;
        end
    end
    
    % faces
    for j = 1:Ny-1
        for i = 1:Nx-1
            v1_idx = Patch_GridVertIdx(i,j);
            v2_idx = Patch_GridVertIdx(i+1,j);
            v3_idx = Patch_GridVertIdx(i+1,j+1);
            v4_idx = Patch_GridVertIdx(i,j+1);
            v5_idx = Patch_CenterVertIdx(i,j);
            
            Faces(face_idx+1,:) = [v5_idx, v1_idx, v2_idx];
            Faces(face_idx+2,:) = [v5_idx, v2_idx, v3_idx];
            Faces(face_idx+3,:) = [v5_idx, v3_idx, v4_idx];
            Faces(face_idx+4,:) = [v5_idx, v4_idx, v1_idx];
            face_idx = face_idx + 4;
        end
    end          
            
end

%% Triangulation -- delete void triangles
solid_idx = find(density >= 0.5); % solid faces
verticesToDelete = setdiff((1:num_vertex)',unique(reshape(Faces(solid_idx,:),[],1)));
verticesToPreserve = setdiff((1:num_vertex)',verticesToDelete);
facesToDelete = setdiff((1:num_tri)',solid_idx);
Vertices(verticesToDelete, :) = [];
Faces(facesToDelete,:) = [];
newVertexIndices = 1:size(Vertices, 1);
vertexMapping = containers.Map(verticesToPreserve, newVertexIndices);
Faces = values(vertexMapping, num2cell(Faces));
Faces = cell2mat(Faces);

%% write obj and stl
fid = fopen(filename_obj,'a');
for i = 1:size(Vertices,1)
    fprintf(fid,'v %f %f %f\n',Vertices(i,1),Vertices(i,2),Vertices(i,3));
end

for i = 1:size(Faces,1)
    fprintf(fid,'f ');
    fprintf(fid,'%d/%d/%d %d/%d/%d %d/%d/%d\n',Faces(i,1),Faces(i,1),Faces(i,1),Faces(i,2),Faces(i,2),Faces(i,2),Faces(i,3),Faces(i,3),Faces(i,3)); 
end

fclose(fid);

stlwrite(triangulation(Faces,Vertices), filename_stl);

end