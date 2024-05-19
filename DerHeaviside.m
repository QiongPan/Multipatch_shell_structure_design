function DerHx = DerHeaviside(x,tau,eta)

%%%%%% The derivative of Heavisede function H(\rho) %%%%%%
%%%%%% H(\rho) = (tanh(tau*eta)+tanh(tau*(\rho-eta)))/(tanh(tau*eta)+tanh(tau*(1-eta)));
% input: 
%   x          - the value of independent variable in H(\rho) 
%   tau        - a parameter controlling the sharpness of the projection function.
%              - tau is gradually doubles to ensure stability.
%   eta        - cut-off value for solid or void.
%              - eta is usually taken as 1/2.
% output: 
%   DerHx      - The derivative of H(\rho)at value x: dH/d{\rho}(\rho=x)

DerHx = tau * (1-tanh(tau*(x-eta))*tanh(tau*(x-eta))) / (tanh(tau*eta) + tanh(tau*(1-eta)));

end