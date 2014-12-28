classdef DataFromOneMolecule < handle
    
    properties (SetAccess = private)
        
        molecule;
        
        cellOfBasisFunctions;
        rhfMBS;
        quambo;
        
    end
    
    methods
        
        function obj = DataFromOneMolecule(mol, basisSetNames)
            import QUAMBO.*;
            
            obj.molecule = mol;
            
            molStr = mol.MoleculeString();
            basisSetMin = basisSetNames.nameMinimalBasisSet;
            basisSetLarge = basisSetNames.LargeBasisSet;
            basisSetLargeAndMin = [basisSetLarge, '_and_', basisSetMin];
            
            % Basis functions
            obj.cellOfBasisFunctions = CellOfBasisFunctions(mol);
            
            % Minimal basis RHF
            matpsi2MBS = MatPsi2(molStr, basisSetMin, 0, 1, [pwd(), '/+QUAMBO']);
            obj.rhfMBS = RHF(RHF.MatPsi2Interface(matpsi2MBS));
            obj.rhfMBS.DoSCF();
            
            % QUAMBO
            obj.quambo = QUAMBO( ...
                QUAMBO.UseMatPsi2(molStr, basisSetLarge, basisSetLargeAndMin));
        end
        
    end
    
end

function cellOfBasisFunctions = CellOfBasisFunctions(mol)
import QUAMBO.BasisFunction;

cellOfBasisFunctions = {};
for iAtom = 1:size(mol.cartesian,1)
    properties.atomicNumber = mol.cartesian(iAtom,1);
    properties.centerXYZ = mol.cartesian(iAtom,2:end);
    if(properties.atomicNumber <= 2) % >=1st row, push 1s
        cellOfBasisFunctions{end+1} = ...
            BasisFunction(Properties1s(properties)); %#ok 1s
    end
    if(properties.atomicNumber <= 10) % >=2nd row, push 2s 2px 2py 2pz
        cellOfBasisFunctions{end+1} = ...
            BasisFunction(Properties2s(properties)); %#ok 2s
        cellOfBasisFunctions{end+1} = ...
            BasisFunction(Properties2px(properties)); %#ok 2px
        cellOfBasisFunctions{end+1} = ...
            BasisFunction(Properties2py(properties)); %#ok 2py
        cellOfBasisFunctions{end+1} = ...
            BasisFunction(Properties2pz(properties)); %#ok 2pz
    end
    if(properties.atomicNumber <= 18) % >=3rd row, push 3s 3px 3py 3pz
        cellOfBasisFunctions{end+1} = ...
            BasisFunction(Properties3s(properties)); %#ok 3s
        cellOfBasisFunctions{end+1} = ...
            BasisFunction(Properties3px(properties)); %#ok 3px
        cellOfBasisFunctions{end+1} = ...
            BasisFunction(Properties3py(properties)); %#ok 3py
        cellOfBasisFunctions{end+1} = ...
            BasisFunction(Properties3pz(properties)); %#ok 3pz
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
