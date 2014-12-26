classdef RHF < handle
    
    properties (SetAccess = private)
        
        overlapMat;
        kineticMat;
        corePotentialMat;
        twoElecIntegrals;
        nuclearRepulsionEnergy;
        numElectrons;
        
        hfEnergy;
        orbital;
        orbitalEnergies;
        
    end
    
    properties (Access = private)
        
        maxSCFIter = 500;
        RMSDensityThreshold = 1e-8;
        MaxDensityThreshold = 1e-6;
        EnergyThreshold = 1e-6;
        
    end
    
    methods
        
        function obj = RHF(properties)
            obj.overlapMat = properties.overlapMat;
            obj.kineticMat = properties.kineticMat;
            obj.corePotentialMat = properties.corePotentialMat;
            obj.twoElecIntegrals = properties.twoElecIntegrals;
            obj.nuclearRepulsionEnergy = properties.nuclearRepulsionEnergy;
            obj.numElectrons = properties.numElectrons;
        end
        
    end
    
    methods (Static)
                
        function properties = MatPsi2Interface(matpsi2)
            properties.overlapMat = matpsi2.Integrals_Overlap();
            properties.kineticMat = matpsi2.Integrals_Kinetic();
            properties.corePotentialMat = matpsi2.Integrals_Potential();
            properties.twoElecIntegrals = matpsi2.Integrals_AllTEIs();
            properties.nuclearRepulsionEnergy = matpsi2.Molecule_NuclearRepulsionEnergy();
            properties.numElectrons = matpsi2.Molecule_NumElectrons();
        end
        
    end
    
end