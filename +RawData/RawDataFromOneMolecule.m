classdef RawDataFromOneMolecule < handle
    
    properties (SetAccess = private)
        
        molecule;
        
        basisFunctions;
        rhfMBS;
        quambo;
        
    end
    
    methods
        
        function obj = RawDataFromOneMolecule(mol, basisSetNames)
            packagePath = what('RawData');
            packagePath = packagePath.path;
            import RawData.*;
            
            obj.molecule = mol;
            
            molStr = mol.MoleculeString();
            basisSetMin = basisSetNames.minimalBasisSet;
            basisSetLarge = basisSetNames.largeBasisSet;
            basisSetLargeAndMin = [basisSetLarge, '_and_', basisSetMin];
            
            % Basis functions
            obj.basisFunctions = CellOfBasisFunctions(mol);
            
            % Minimal basis RHF
            matpsi2MBS = MatPsi2(molStr, basisSetMin, 0,1,packagePath);
            obj.rhfMBS = RHF(RHF.MatPsi2Interface(matpsi2MBS));
            obj.rhfMBS.DoSCF();
            
            % QUAMBO
            obj.quambo = QUAMBO(QUAMBO.UseMatPsi2( ...
                molStr, basisSetLarge, basisSetLargeAndMin));
        end
        
    end
    
end

function basisFunctions = CellOfBasisFunctions(mol)
import RawData.BasisFunction;

basisFunctions = {};
for iAtom = 1:size(mol.cartesian,1)
    properties.atomicNumber = mol.cartesian(iAtom,1);
    properties.centerXYZ = mol.cartesian(iAtom,2:end);
    if(properties.atomicNumber <= 2) % >=1st row, push 1s
        basisFunctions{end+1} = ...
            BasisFunction(Properties1s(properties)); %#ok 1s
    end
    if(properties.atomicNumber <= 10) % >=2nd row, push 2s 2px 2py 2pz
        basisFunctions{end+1} = ...
            BasisFunction(Properties2s(properties)); %#ok 2s
        basisFunctions{end+1} = ...
            BasisFunction(Properties2px(properties)); %#ok 2px
        basisFunctions{end+1} = ...
            BasisFunction(Properties2py(properties)); %#ok 2py
        basisFunctions{end+1} = ...
            BasisFunction(Properties2pz(properties)); %#ok 2pz
    end
    if(properties.atomicNumber <= 18) % >=3rd row, push 3s 3px 3py 3pz
        basisFunctions{end+1} = ...
            BasisFunction(Properties3s(properties)); %#ok 3s
        basisFunctions{end+1} = ...
            BasisFunction(Properties3px(properties)); %#ok 3px
        basisFunctions{end+1} = ...
            BasisFunction(Properties3py(properties)); %#ok 3py
        basisFunctions{end+1} = ...
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
