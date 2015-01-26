classdef MLHueckel < handle
    
    properties (SetAccess = private)
        
        numElectrons;
        nucRepEnergy;
        
        overlapMat;
        coreHamilt;
        fockMat;
        
    end
    
    methods
        
        function obj = MLHueckel(input)
            obj.numElectrons = input.numElectrons;
            obj.nucRepEnergy = input.nucRepEnergy;
            obj.overlapMat = input.overlapMat;
            obj.coreHamilt = input.coreHamilt;
            obj.fockMat = input.fockMat;
        end
        
        function solution = Solve(obj)
            invSHalf = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);
            fockMatOrtho = invSHalf * obj.fockMat * invSHalf;
            [eigVectors, eigValues] = eig(fockMatOrtho);
            [solution.orbitalEnergies, order] = sort(diag(eigValues));
            solution.orbital = invSHalf * eigVectors(:, order);
            occOrb = solution.orbital(:,1:obj.numElectrons/2);
            solution.densKernelMat = occOrb * occOrb';
            solution.totalEnergy = reshape(solution.densKernelMat, 1,[]) ...
                * reshape(obj.coreHamilt+obj.fockMat, [],1) ...
                + obj.nucRepEnergy;
        end
        
    end
    
end
