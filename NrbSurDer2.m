/**
 * @file NrbSurDer2.m
 * @author Qiong Pan (PQ2019@mail.ustc.edu.cn)
 * @brief  Non-Zero NURBS basis functions and their first, second partial derivatives in a surface point
 * Orginal work written by Qiong Pan in Matlab.
 * @date 2021-11-1
 */

function [IdxU,IdxV, N, Ns, Nt, Nst, Nss, Ntt] = NrbSurDer2(m,n,p,q,u,v,U,V,Weights)

%%%%%% The 0-th. first and second partial derivatives of (non-zero) NURBS basis functions at (u,v) %%%%%% 
% input: 
%   m          - number of control points - 1 (U-direction)/(s-direction)
%   n          - number of control points - 1 (V-direction)/(t-direction)
%   p          - spline degree   (U-direction)
%   q          - spline degree   (V-direction)
%   u          - parameter value (U-direction)
%   v          - parameter value (V-direction)
%   U          - knot vector     (U-direction)
%   V          - knot vector     (V-direction)
%   Weights    - weights for NURBS basis functions
%              - dim = (m+1) * (n+1)
% output: 
%   IdxU       - Indices of Non-zero basis functions at u in u-direction, a vector 
%              - dim = (p+1) * 1
%   IdxV       - Indices of Non-zero basis functions at v in v-direction, a vector 
%              - dim = (q+1) * 1
%   N          - Non-zero basis functions in (u,v), a matrix
%              - dim = (p+1) * (q+1)
%   N_s        - first partial derivative w.r.t. u for Non-zero basis functions in (u,v), a matrix
%              - dim = (p+1) * (q+1)
%   N_t        - first partial derivative w.r.t. v for Non-zero basis functions in (u,v), a matrix
%              - dim = (p+1) * (q+1)
%   N_st       - second partial derivative w.r.t. u and v for Non-zero basis functions in (u,v), a matrix
%              - dim = (p+1) * (q+1)
%   N_ss       - second partial derivative w.r.t. u for Non-zero basis functions in (u,v), a matrix
%              - dim = (p+1) * (q+1)
%   N_tt       - second partial derivative w.r.t. v for Non-zero basis functions in (u,v), a matrix
%              - dim = (p+1) * (q+1)

r1 = FindSpan(m, p, u, U);
r2 = FindSpan(n, q, v, V);

IdxU = NonZeroBasisIdx(m, p, u, U);
IdxV = NonZeroBasisIdx(n, q, v, V);   

Ders2   = BasisFuncDer(r1, p, u, U, 2);    % dim = 3 * (p+1)
Dert2   = BasisFuncDer(r2, q, v, V, 2);    % dim = 3 * (q+1)  

WIJ     = Weights(IdxU+1,IdxV+1);
WNIJ     = Ders2(1,:)' * Dert2(1,:) .* WIJ;    	
WNIJ_s   = Ders2(2,:)' * Dert2(1,:) .* WIJ;   	
WNIJ_t   = Ders2(1,:)' * Dert2(2,:) .* WIJ;    	   	
WNIJ_ss  = Ders2(3,:)' * Dert2(1,:) .* WIJ; 	
WNIJ_tt  = Ders2(1,:)' * Dert2(3,:) .* WIJ; 
WNIJ_st  = Ders2(2,:)' * Dert2(2,:) .* WIJ; 

W = sum(WNIJ, 'all');
W_s = sum(WNIJ_s, 'all');
W_t = sum(WNIJ_t, 'all');
W_ss = sum(WNIJ_ss, 'all');
W_tt = sum(WNIJ_tt, 'all');
W_st = sum(WNIJ_st, 'all');

N   = WNIJ / W;                                         	% rational basis of S(u,v) 
Ns  = (WNIJ_s  - W_s  * N) / W;                         	% Ss = (As-S*Ws)/W
Nt  = (WNIJ_t  - W_t  * N) / W;                             % St = (At-S*Wt)/W
Nss = (WNIJ_ss - W_ss * N - 2 * W_s * Ns) / W;          	% Sss = (Ass-Wss*S-2*Ws*Ss)/W
Ntt = (WNIJ_tt - W_tt * N - 2 * W_t * Nt) / W;          	% Stt = (Att-Wtt*S-2*Wt*St)/W
Nst = (WNIJ_st - W_st * N - W_s * Nt - W_t * Ns) / W;       % Sst = (Ast-Wst*S-Ws*St-Wt*Ss)/W

end
