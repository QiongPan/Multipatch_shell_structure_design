function [c1,c2,c3,c4] = ComputeCircumference(p,q,U,V,HCtrlpts)

%%%%%% Compute the length of four edges of the NRBS surface %%%%%%
% input: 
%   p          - spline degree   (U-direction)
%   q          - spline degree   (V-direction)
%   U          - knot vector     (U-direction)
%   V          - knot vector     (V-direction)
%   Ctrlpts    - control points (m+1) * (n+1) * 4
% output: 
%   c1         - length of u-curve  S(u,v=0)
%   c2         - length of v-curve  S(u=1,v)
%   c3         - length of u-curve  S(u,v=1)
%   c4         - length of v-curve  S(u-0,v)

global ng_u ng_v;
ne_s = 100;
ne_t = 100;

Length_All  = zeros(ne_s, 2);
Width_All   = zeros(ne_t, 2);
UEle  	= linspace(U(1),U(end),ne_s+1);
VEle   	= linspace(V(1),V(end),ne_t+1);
uins 	= setdiff(UEle, U);
vins  	= setdiff(VEle, V);
[UAnalysis, VAnalysis, HCtrlpts] = KnotInsert(uins,vins,U,V,HCtrlpts); 
Ctrlpts             = HCtrlpts;
Ctrlpts(:,:,1:3)  	= Ctrlpts(:,:,1:3)./Ctrlpts(:,:,4);
[mAnalysis,nAnalysis,~]    = size(Ctrlpts);
mAnalysis                  = mAnalysis - 1; 
nAnalysis                  = nAnalysis - 1;
Weights = reshape(Ctrlpts(:,:,4),mAnalysis+1,nAnalysis+1);

for i = 1:ne_s
    [svec, w1] = GaussInt(UEle(i),UEle(i+1),ng_u);
    tvec = [VEle(1); VEle(end)];
    for j = 1:2    
        v = tvec(j);
        for k = 1:numel(w1) 
            u = svec(k);
            [IdxU, IdxV, ~, Ns, ~, ~, ~, ~] = NrbSurDer2(mAnalysis,nAnalysis,p,q,u,v,UAnalysis,VAnalysis,Weights);
            PIJ = Ctrlpts(IdxU+1,IdxV+1,1:3);
            Ss = reshape(sum(PIJ.*Ns, [1,2]),1,3);
            Length_All(k,j) = Length_All(k,j) + w1(k) * norm(Ss); % guass integration
        end   
    end
end

for j = 1:ne_t
    [tvec, w2] = GaussInt(VEle(j),VEle(j+1),ng_v);
    svec = [UEle(1); UEle(end)];
    for i = 1:2    
        u = svec(i);
        for k = 1:numel(w2) 
            v = tvec(k);
            [IdxU, IdxV, ~, ~, Nt, ~, ~, ~] = NrbSurDer2(mAnalysis,nAnalysis,p,q,u,v,UAnalysis,VAnalysis,Weights);
            PIJ = Ctrlpts(IdxU+1,IdxV+1,1:3);
            St = reshape(sum(PIJ.*Nt, [1,2]),1,3);
            Width_All(k,i) = Width_All(k,i) + w2(k) * norm(St);
        end    
    end
end

c1 = sum(Length_All(:,1));
c2 = sum(Width_All(:,2));
c3 = sum(Length_All(:,2));
c4 = sum(Width_All(:,1));
end