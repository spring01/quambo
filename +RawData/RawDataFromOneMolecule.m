classdef RawDataFromOneMolecule < handle
    
    properties (SetAccess = private)
        
        molecule;
        
        basisFunctions;
        rhfAMBO;
        rhfAO;
        
        quambo;
        
    end
    
    methods
        
        function obj = RawDataFromOneMolecule(molecule, basisSetNames)
            packagePath = what('RawData');
            packagePath = packagePath.path;
            import RawData.*;
            
            obj.molecule = molecule;
            
            molStr = molecule.MoleculeString();
            basisSetAO = basisSetNames.basisSetAO;
            basisSetAMBO = basisSetNames.basisSetAMBO;
            auxiliaryBasisSet = basisSetNames.auxiliaryBasisSet;
            
            % Basis functions
            obj.basisFunctions = CellOfBasisFunctions(molecule);
            
            % Minimal basis RHF
            matpsi2AMBO = MatPsi2(molStr, basisSetAMBO, 0,1,packagePath);
            matpsi2AMBO.JK_Initialize('DFJK', auxiliaryBasisSet);
            obj.rhfAMBO = RHF(RHF.MatPsi2Interface(matpsi2AMBO));
            matpsi2AO = MatPsi2(molStr, basisSetAO);
            matpsi2AO.JK_Initialize('DFJK', auxiliaryBasisSet);
            obj.rhfAO = RHF(RHF.MatPsi2Interface(matpsi2AO));
            
            % QUAMBO
            obj.quambo = QUAMBO(QUAMBO.UseMatPsi2(molStr, basisSetNames));
        end
        
    end
    
end

function basisFunctions = CellOfBasisFunctions(mol)
import RawData.BasisFunction;

basisFunctions = {};
for iAtom = 1:size(mol.cartesian,1)
    properties.atomicNumber = mol.cartesian(iAtom,1);
    properties.centerXYZ = mol.cartesian(iAtom,2:end);
    basisFunctions{end+1} = ... % always push 1s
        BasisFunction(Properties1s(properties)); %#ok
    if(properties.atomicNumber > 2) % >=2nd row, push 2s 2px 2py 2pz
        basisFunctions{end+1} = ...
            BasisFunction(Properties2s(properties)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Properties2px(properties)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Properties2py(properties)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Properties2pz(properties)); %#ok
    end
    if(properties.atomicNumber > 10) % >=3rd row, push 3s 3px 3py 3pz
        basisFunctions{end+1} = ...
            BasisFunction(Properties3s(properties)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Properties3px(properties)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Properties3py(properties)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Properties3pz(properties)); %#ok
    end
    if(properties.atomicNumber > 18)
        throw(MException('BasisSet:BasisSet','Atomic number too large.'));
    end
end
end

function properties = Properties1s(properties)
properties.principleQNum = 1;
properties.angularQNum = 0;
properties.cartesianQNum = 0;
end

function properties = Properties2s(properties)
properties.principleQNum = 2;
properties.angularQNum = 0;
properties.cartesianQNum = 0;
end

function properties = Properties2px(properties)
properties.principleQNum = 2;
properties.angularQNum = 1;
properties.cartesianQNum = 1;
end

function properties = Properties2py(properties)
properties.principleQNum = 2;
properties.angularQNum = 1;
properties.cartesianQNum = 2;
end

function properties = Properties2pz(properties)
properties.principleQNum = 2;
properties.angularQNum = 1;
properties.cartesianQNum = 3;
end

function properties = Properties3s(properties)
properties.principleQNum = 3;
properties.angularQNum = 0;
properties.cartesianQNum = 0;
end

function properties = Properties3px(properties)
properties.principleQNum = 3;
properties.angularQNum = 1;
properties.cartesianQNum = 1;
end

function properties = Properties3py(properties)
properties.principleQNum = 3;
properties.angularQNum = 1;
properties.cartesianQNum = 2;
end

function properties = Properties3pz(properties)
properties.principleQNum = 3;
properties.angularQNum = 1;
properties.cartesianQNum = 3;
end
