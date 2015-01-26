classdef RawData1Mol < handle
    
    properties (SetAccess = private)
        
        molecule;
        
        basisFunctions;
        
        mlHueckelTarget;
        
    end
    
    methods
        
        function obj = RawData1Mol(molecule, basisSetNames)
            import RawData.*;
            
            obj.molecule = molecule;
            
            obj.basisFunctions = CellOfBasisFunctions(molecule);
            
            molStr = molecule.MoleculeString();
            packagePath = what('RawData');
            packagePath = packagePath.path;
            basisSetAO = basisSetNames.basisSetAO;
            basisSetAOandAMBO = [basisSetNames.basisSetAO, '_and_', basisSetNames.basisSetAMBO];
            matpsi2AO = MatPsi2(molStr, basisSetAO);
            matpsi2AO.RHF_DoSCF();
            matpsi2AOandAMBO = MatPsi2(molStr, basisSetAOandAMBO, 0, 1, packagePath);
            
            % QUAMBO
            quambo = QUAMBO(QUAMBO.UseMatPsi2(matpsi2AO, matpsi2AOandAMBO));
            ao2quambo = quambo.AOtoQUAMBO;
            input.numElectrons = matpsi2AO.Molecule_NumElectrons();
            input.nucRepEnergy = matpsi2AO.Molecule_NuclearRepulsionEnergy();
            input.overlapMat = ao2quambo' * matpsi2AO.Integrals_Overlap() * ao2quambo;
            input.coreHamilt = ao2quambo' * matpsi2AO.RHF_CoreHamiltonian() * ao2quambo;
            input.fockMat = ao2quambo' * matpsi2AO.RHF_FockMatrix() * ao2quambo;
            obj.mlHueckelTarget = MLHueckel(input);
        end
        
    end
    
end


function basisFunctions = CellOfBasisFunctions(mol)
import RawData.BasisFunction;

basisFunctions = {};
for iAtom = 1:size(mol.cartesian,1)
    input.atomicNumber = mol.cartesian(iAtom,1);
    input.centerXYZ = mol.cartesian(iAtom,2:end);
    basisFunctions{end+1} = ... % always push 1s
        BasisFunction(Input1s(input)); %#ok
    if(input.atomicNumber > 2) % >=2nd row, push 2s 2px 2py 2pz
        basisFunctions{end+1} = ...
            BasisFunction(Input2s(input)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Input2px(input)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Input2py(input)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Input2pz(input)); %#ok
    end
    if(input.atomicNumber > 10) % >=3rd row, push 3s 3px 3py 3pz
        basisFunctions{end+1} = ...
            BasisFunction(Input3s(input)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Input3px(input)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Input3py(input)); %#ok
        basisFunctions{end+1} = ...
            BasisFunction(Input3pz(input)); %#ok
    end
    if(input.atomicNumber > 18)
        throw(MException('BasisSet:BasisSet','Atomic number too large.'));
    end
end
end

function input = Input1s(input)
input.principalQNum = 1;
input.angularQNum = 0;
input.cartesianQNum = 0;
end

function input = Input2s(input)
input.principalQNum = 2;
input.angularQNum = 0;
input.cartesianQNum = 0;
end

function input = Input2px(input)
input.principalQNum = 2;
input.angularQNum = 1;
input.cartesianQNum = 1;
end

function input = Input2py(input)
input.principalQNum = 2;
input.angularQNum = 1;
input.cartesianQNum = 2;
end

function input = Input2pz(input)
input.principalQNum = 2;
input.angularQNum = 1;
input.cartesianQNum = 3;
end

function input = Input3s(input)
input.principalQNum = 3;
input.angularQNum = 0;
input.cartesianQNum = 0;
end

function input = Input3px(input)
input.principalQNum = 3;
input.angularQNum = 1;
input.cartesianQNum = 1;
end

function input = Input3py(input)
input.principalQNum = 3;
input.angularQNum = 1;
input.cartesianQNum = 2;
end

function input = Input3pz(input)
input.principalQNum = 3;
input.angularQNum = 1;
input.cartesianQNum = 3;
end