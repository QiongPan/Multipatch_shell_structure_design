function Hx = Heaviside(x,tau,eta)

%%%%%% The Heavisede function H(\rho) %%%%%%
% input: 
%   x          - the value of independent variable in H(\rho) 
%   tau        - a parameter controlling the sharpness of the projection function.
%              - tau is gradually doubles to ensure stability.
%   eta        - cut-off value for solid or void.
%              - eta is usually taken as 1/2.
% output: 
%   Hx         - the values of the Heaviside function at x: H(\rho=x)

Hx = (tanh(tau*eta)+tanh(tau*(x-eta)))/(tanh(tau*eta)+tanh(tau*(1-eta)));

end
