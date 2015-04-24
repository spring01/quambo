classdef ModifiedRHF < RHF
    
    properties (SetAccess = private)
        
        overlapMatModifier;
        kineticMatModifier;
        corePotentialMatsModifier;
        coulombMatModifier;
        exchangeMatModifier;
        
    end
    
    methods
        
        function obj = ModifiedRHF(matpsi2, rhfQuambo)
            rhfll = RHF.MatPsi2Interface(matpsi2);
            rhfll.SCF();
            obj = obj@RHF(rhfll);
            
            densVec = rhfQuambo.densvec;
            teiForCoulomb = reshape(rhfll.twoElecIntegrals, length(densVec), []);
            teiForExchange = reshape(permute(rhfll.twoElecIntegrals, [1 3 2 4]), ...
                length(densVec), []);
            llcoulombMat = reshape(teiForCoulomb * densVec, size(obj.overlapMat));
            llexchangeMat = reshape(teiForExchange * densVec, size(obj.overlapMat));
            
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
            
            obj.coulombMatModifier = rhfQuambo.coulombMat ./ llcoulombMat;
            obj.coulombMatModifier(isnan(obj.coulombMatModifier)) = 1;
            obj.coulombMatModifier(abs(obj.coulombMatModifier) > thres) = 1;
            
            obj.exchangeMatModifier = rhfQuambo.exchangeMat ./ llexchangeMat;
            obj.exchangeMatModifier(isnan(obj.exchangeMatModifier)) = 1;
            obj.exchangeMatModifier(abs(obj.exchangeMatModifier) > thres) = 1;
        end
        
        [hfEnergy, orbital, orbitalEnergies, iter] = SCF(obj);
        
    end
    
end