classdef FeatureLabelCoreHamiltOffDiag < FeatureLabel.FeatureLabelBase
    
    methods
        
        function obj = FeatureLabelCoreHamiltOffDiag(rawData1Mol, ind1, ind2)
            obj.identifier = obj.OffDiagIdentifier(rawData1Mol, ind1, ind2);
            
            % feature
            obj.feature = obj.FeatureCoreHamiltOffDiag(rawData1Mol, ind1, ind2);
            
            % label
            obj.label = rawData1Mol.mltb.coreHamilt(ind1,ind2);
        end
        
    end
    
    methods (Access = private)
        
        function feature = FeatureCoreHamiltOffDiag(obj, rawData1Mol, ind1, ind2)
            func1 = rawData1Mol.basisFunctions{ind1};
            func2 = rawData1Mol.basisFunctions{ind2};
            feature = [ ...
                func1.atomicNumber, ...
                func1.principalQNum, ...
                func1.angularQNum, ...
                func1.cartesianQNum, ...
                func2.atomicNumber, ...
                func2.principalQNum, ...
                func2.angularQNum, ...
                func2.cartesianQNum, ...
                func2.centerXYZ - func1.centerXYZ];
        end
        
    end
    
end
