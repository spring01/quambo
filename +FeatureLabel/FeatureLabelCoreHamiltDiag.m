classdef FeatureLabelCoreHamiltDiag < FeatureLabel.FeatureLabelBase
    
    methods
        
        function obj = FeatureLabelCoreHamiltDiag(rawData1Mol, ind)
            obj.identifier = obj.DiagIdentifier(rawData1Mol, ind);
            
            % feature
            obj.feature = obj.FeatureCoreHamiltDiag(rawData1Mol, ind);
            
            % label
            obj.label = rawData1Mol.mltb.coreHamilt(ind,ind);
        end
        
    end
    
    methods (Access = private)
        
        function feature = FeatureCoreHamiltDiag(obj, rawData1Mol, ind)
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
