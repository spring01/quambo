classdef BasisFunction < handle
    
    properties (SetAccess = private)
        
        atomicNumber; % Z
        centerXYZ; % [x y z]
        
        principalQNum; % n
        angularQNum; % l
        cartesianQNum; % lc: cartesian version of lz; x->1, y->2, z->3
        
    end
    
    methods
        
        function obj = BasisFunction(input)
            obj.atomicNumber = input.atomicNumber;
            obj.centerXYZ = input.centerXYZ;
            obj.principalQNum = input.principalQNum;
            obj.angularQNum = input.angularQNum;
            obj.cartesianQNum = input.cartesianQNum;
        end
        
    end
    
end
