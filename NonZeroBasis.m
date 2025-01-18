/**
 * @file KnotInsert.m
 * @author Qiong Pan (PQ2019@mail.ustc.edu.cn)
 * @brief  Non-Zero NURBS basis finctions for a NURBS surface
 * Orginal work written by Qiong Pan in Matlab.
 * @date 2021-11-1
 */

function [Val, IdxU, IdxV] = NonZeroBasis(m, n, p, q, u, v, U, V, Weights)

%%%%%% Non-zero NURBS basis functions and indices for a NURBS surface in (u,v) %%%%%%
% input: 
%   m          - number of control points - 1 (U-direction)
%   n          - number of control points - 1 (V-direction)
%   p          - spline degree   (U-direction)
%   q          - spline degree   (V-direction)
%   u          - parameter value (U-direction)
%   v          - parameter value (V-direction)
%   U          - uknot vector    (U-direction)
%   V          - vknot vector    (V-direction)
%   Weights    - weights for NURBS basis functions
%              - dim = (m+1) * (n+1)
% output: 
%   Val        - Non-zero basis functions in (u,v), a matrix
%              - dim = (p+1) * (q+1)
%   IdxU       - Indices of Non-zero basis functions at u in u-direction, a vector 
%              - dim = (p+1) * 1
%   IdxV       - Indices of Non-zero basis functions at v in v-direction, a vector 
%              - dim = (q+1) * 1

r1 = FindSpan(m, p, u, U); 
r2 = FindSpan(n, q, v, V);

IdxU = NonZeroBasisIdx(m, p, u, U);
IdxV = NonZeroBasisIdx(n, q, v, V);

Bu = BasisFunc(r1, p, u, U);
Bv = BasisFunc(r2, q, v, V);

NIJ = Bu * Bv';
WIJ = Weights(IdxU+1,IdxV+1);
NIJ = NIJ.*WIJ/sum(NIJ.*WIJ,'all');
Val = NIJ;
end
