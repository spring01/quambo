classdef BasisFunction < handle
    
    properties (SetAccess = private)
        
        atomicNumber; % Z
        centerXYZ; % [x y z]
        
        principleQNum; % n
        angularQNum; % l
        cartesianQNum; % lc: cartesian version of lz; x->1, y->2, z->3
        
    end
    
    methods
        
        function obj = BasisFunction(properties)
            obj.atomicNumber = properties.atomicNumber;
            obj.centerXYZ = properties.centerXYZ;
            obj.principleQNum = properties.principleQNum;
            obj.angularQNum = properties.angularQNum;
            obj.cartesianQNum = properties.cartesianQNum;
        end
        
    end
    
end
