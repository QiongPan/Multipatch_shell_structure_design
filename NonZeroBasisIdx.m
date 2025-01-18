/**
 * @file NonZeroBasisIdx.m
 * @author Qiong Pan (PQ2019@mail.ustc.edu.cn)
 * @brief Indices of Non-Zero NURBS basis functions for a NURBS surface
 * Orginal work written by Qiong Pan in Matlab.
 * @date 2021-11-1
 */

function idx = NonZeroBasisIdx(n, p, u, U)

%%%%%% Numbering of basis functions for NURBS %%%%%%
% input:
%   n          - number of control points - 1
%   p          - spline degree
%   u          - parameter
%   U          - knot vector
% output:
%   idx        - Indices of the basis functions that are nonvanishing at u, a vector
%              - dim = (p+1) * 1

r = FindSpan(n, p, u, U);   % knot span index 
idx = r - p + (0:p)';       % idx = [r-p, r-p+1, ..., r], start from 0
  
end
