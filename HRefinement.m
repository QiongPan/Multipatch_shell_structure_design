/**
 * @file HRefinement.m
 * @author Qiong Pan (PQ2019@mail.ustc.edu.cn)
 * @brief  compute the relation matrix between corse and fine control meshes for a NURBS surface
 * Orginal work written by Qiong Pan in Matlab.
 * @date 2021-11-1
 */

function [UNew,VNew,UTransMat,VTransMat] = HRefinement(p,q,m,n,U,V,Uins,Vins)

%%%%%% Compute the transformation matrix between control points of the former and refined surfaces  %%%%%%         
% input:
%   p          - the degree of basis functions for u-curves
%   q          - the degree of basis functions for v-curves
%   m          - number of control points - 1, in u-direction, starting from 0
%   n          - number of control points - 1, in v-direction, starting from 0
%   U          - uknot vector
%   V          - vknot vector
%   Uins       - knot sequence to be inserted into U in  u-direction, a vector
%   Vins       - knot sequence to be inserted into V in  v-direction, a vector
% output:
%   UNew       - uknot vector of the new surface
%   VNew       - vknot vector of the new surface
%   UTransMat  - Linear transformation matrix between control points of the former and refined u-curves
%   VTransMat  - Linear transformation matrix between control points of the former and refined v-curves

UTransMat = eye(m+1);
VTransMat = eye(n+1);
UNew = U;
VNew = V;
m_ = m;
n_ = n;
for i = 1:numel(Uins)
    ubar = Uins(i);
    r = FindSpan(m_, p, ubar, UNew);  
    UNew = sort([UNew, ubar]);
    TransMat = zeros(m_+1,m_+2);
    TransMat(1:r-p,1:r-p) = eye(r-p);
    TransMat(r+2:end,r+3:end) = eye(m_-r);
    for ii = (r-p+1):(r+1)
        TransMat(ii,ii) = (ubar-UNew(ii)) / (UNew(ii+p+1)-UNew(ii));
        TransMat(ii,ii+1) = (UNew(ii+p+2)-ubar) / (UNew(ii+p+2)-UNew(ii+1));
    end
    UTransMat = UTransMat * TransMat;
    m_ = m_ + 1;
end

for j = 1:numel(Vins)
    vbar = Vins(j);
    r = FindSpan(n_, q, vbar, VNew);  
    VNew = sort([VNew, vbar]);
    TransMat = zeros(n_+1,n_+2);
    TransMat(1:r-q,1:r-q) = eye(r-q);
    TransMat(r+2:end,r+3:end) = eye(n_-r);
    for jj = (r-q+1):(r+1)
        TransMat(jj,jj) = (vbar-VNew(jj)) / (VNew(jj+q+1)-VNew(jj));
        TransMat(jj,jj+1) = (VNew(jj+q+2)-vbar) / (VNew(jj+q+2)-VNew(jj+1));
    end
    VTransMat = VTransMat * TransMat;
    n_ = n_ + 1;
end


end
