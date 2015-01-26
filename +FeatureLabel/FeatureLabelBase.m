classdef (Abstract) FeatureLabelBase < handle
    
    properties (SetAccess = protected)
        
        identifier;
        
        feature;
        label;
        
    end
    
    methods (Access = protected)
        
        function identifier = DiagIdentifier(~, rawData1Mol, ind)
            func = rawData1Mol.basisFunctions{ind};
            identifier ...
                = num2str( ...
                1000000 .* func.atomicNumber + ...
                10000 .* func.principalQNum + ...
                100 .* func.angularQNum + ...
                1 .* func.cartesianQNum);
        end
        
        function identifier = OffDiagIdentifier(obj, rawData1Mol, ind1, ind2)
            id1 = obj.DiagIdentifier(rawData1Mol, ind1);
            id2 = obj.DiagIdentifier(rawData1Mol, ind2);
            if(str2double(id1) < str2double(id2))
                temp = id1;
                id1 = id2;
                id2 = temp;
            end
            identifier = [id1, id2];
        end
        
    end
    
end
