classdef FeatureLabelOverlap < QUAMBO.FeatureLabelBase
    
    methods
        
        function obj = FeatureLabelOverlap(data1mol, i, j)
            feature1 = data1mol.rhfMBS.overlapMat(i,j);
            feature2 = obj.FeatureOverlap( ...
                data1mol.cellOfBasisFunctions{i}, ...
                data1mol.cellOfBasisFunctions{j});
            
            obj.feature = [feature1, feature2];
            obj.label = data1mol.quambo.overlapQUAMBO(i, j);
        end
        
    end
    
    methods (Access = private)
        
        function feature = FeatureOverlap(obj, bfi, bfj)
        end
        
    end
    
end
