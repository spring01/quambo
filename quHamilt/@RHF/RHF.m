classdef RHF < handle
    
    properties (SetAccess = protected)
        
        overlapMat;
        kineticMat;
        corePotentialMats;
        twoElecIntegrals;
        nucRepEnergy;
        numElectrons;
        
        hfEnergy;
        orbital;
        orbitalEnergies;
        densvec
        
        coulombMat;
        exchangeMat;
        
    end
    
    properties (Access = protected)
        
        maxSCFIter = 500;
        RMSDensityThreshold = 1e-8;
        MaxDensityThreshold = 1e-6;
        EnergyThreshold = 1e-6;
        
    end
    
    methods
        
        function obj = RHF(properties)
            obj.overlapMat = properties.overlapMat;
            obj.kineticMat = properties.kineticMat;
            obj.corePotentialMats = properties.corePotentialMats;
            obj.twoElecIntegrals = properties.twoElecIntegrals;
            obj.nucRepEnergy = properties.nucRepEnergy;
            obj.numElectrons = properties.numElectrons;
        end
        
        [hfEnergy, orbital, orbitalEnergies, iter] = SCF(obj);
        
    end
    
    methods (Static)
                
        function rhf = MatPsi2Interface(matpsi2)
            properties.overlapMat = matpsi2.Integrals_Overlap();
            properties.kineticMat = matpsi2.Integrals_Kinetic();
            properties.corePotentialMats = matpsi2.Integrals_PotentialEachCore();
            properties.twoElecIntegrals = matpsi2.Integrals_AllTEIs();
            properties.nucRepEnergy = matpsi2.Molecule_NucRepEnergy();
            properties.numElectrons = matpsi2.Molecule_NumElectrons();
            rhf = RHF(properties);
        end
        
    end
    
end