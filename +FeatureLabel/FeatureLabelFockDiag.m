classdef FeatureLabelFockDiag < FeatureLabel.FeatureLabelBase
    
    methods
        
        function obj = FeatureLabelFockDiag(rawData1Mol, ind)
            obj.identifier = obj.DiagIdentifier(rawData1Mol, ind);
            
            % feature
            obj.feature = obj.FeatureFockDiag(rawData1Mol, ind);
            
            % label
            obj.label = rawData1Mol.mltb.fockMat(ind,ind);
        end
        
    end
    
    methods (Access = private)
        
        function feature = FeatureFockDiag(obj, rawData1Mol, ind)
            func = rawData1Mol.basisFunctions{ind};
            feature = [ ...
                func.atomicNumber, ...
                func.principalQNum, ...
                func.angularQNum, ...
                func.cartesianQNum, ...
                func.centerXYZ];
        end
        
    end
    
end
