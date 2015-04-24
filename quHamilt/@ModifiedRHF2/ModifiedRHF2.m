classdef ModifiedRHF2 < RHF
    
    properties (SetAccess = private)
        
        overlapMatModifier;
        kineticMatModifier;
        corePotentialMatsModifier;
        teiModifier;
        
    end
    
    methods
        
        function obj = ModifiedRHF2(matpsi2, rhfQuambo)
            rhfll = RHF.MatPsi2Interface(matpsi2);
            rhfll.SCF();
            obj = obj@RHF(rhfll);
            
            thres = 1e3;
            
            obj.overlapMatModifier = rhfQuambo.overlapMat ./ obj.overlapMat;
            obj.overlapMatModifier(isnan(obj.overlapMatModifier)) = 1;
            obj.overlapMatModifier(abs(obj.overlapMatModifier) > thres) = 1;
            
            obj.kineticMatModifier = rhfQuambo.kineticMat ./ obj.kineticMat;
            obj.kineticMatModifier(isnan(obj.kineticMatModifier)) = 1;
            obj.kineticMatModifier(abs(obj.kineticMatModifier) > thres) = 1;
            
            obj.corePotentialMatsModifier = rhfQuambo.corePotentialMats ./ obj.corePotentialMats;
            obj.corePotentialMatsModifier(isnan(obj.corePotentialMatsModifier)) = 1;
            obj.corePotentialMatsModifier(abs(obj.corePotentialMatsModifier) > thres) = 1;
            
            obj.teiModifier = rhfQuambo.twoElecIntegrals ./ rhfll.twoElecIntegrals;
            obj.teiModifier(isnan(obj.teiModifier)) = 1;
            obj.teiModifier(abs(obj.teiModifier) > thres) = 1;
        end
        
        [hfEnergy, orbital, orbitalEnergies, iter] = SCF(obj);
        
    end
    
end