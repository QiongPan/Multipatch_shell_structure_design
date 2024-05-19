function mult = FindMultiplicity(u, U)

%%%%%% Multiplicity of the parameter u in knot vector U %%%%%%
% input: 
%   u          - parameter u \in [u_i, u_{i+1})
%   U          - knot vector  
% output: 
%   mult       - the multiplicity of u in U

    if (u>U(end) || u<U(1))
        error ('the value is outside the knot interval of definition');
    end
    
    mult = 0;
    for i=1:length(U)
        if u > U(i)
            continue;
        elseif(u == U(i))
            mult = mult + 1;        
        else
            break;
        end
    end
    
end