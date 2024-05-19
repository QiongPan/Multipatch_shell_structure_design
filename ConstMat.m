function [H, D] = ConstMat(E, mu)

%%%%%% The constant matrices H and D used to compute the stiffness matrix K %%%%%%
%%%%%% K = \int_{\Omega} B^{T}* D * B d\Omega 
%%%%%% B is the stain-displacement matrix. B = T * H * {\Gamma} * R. 
%%%%%% D is the stain-stress matrix.
% input:
%   E          - Young's modulus
%   nu         - Poisson's ratio
% output: 
%   H          - the transformation matrix between strain vector and displacement derivatives w.r.t. x,y,z
%              - ({\epsilon}_{xx}, {\epsilon}_{yy}, ... {\gamma}_{xz}) = H * (du/dx, du/dy, du/dz, dv/dx, ... dw/dz)
%              - dim = 6*9
%   D          - the stain-stress matrix
%              - dim = 6*6

a = E/(1 - mu^2);
b = mu * a;
c = E / 2 / (1 + mu);
d = 5 / 6 * c;

H = [1 0 0 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0 1;
     0 1 0 1 0 0 0 0 0;
     0 0 0 0 0 1 0 1 0;
     0 0 1 0 0 0 1 0 0];

D = [a b 0 0 0 0;
     b a 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 c 0 0;
     0 0 0 0 d 0;
     0 0 0 0 0 d];
 
end