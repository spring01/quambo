function input = UseMatPsi2(matpsi2AO, matpsi2AOandAMBO)

input.numElectrons = matpsi2AO.Molecule_NumElectrons();
input.AOtoMO = matpsi2AO.RHF_Orbital();

atomicNumbers = matpsi2AO.Molecule_AtomicNumbers();
centerNumFunctionsAO = ...
    CenterNumFunctions(atomicNumbers, matpsi2AO.BasisSet_FunctionToCenter());
centerNumFunctionsAMBO = ...
    CenterNumFunctions(atomicNumbers, matpsi2AOandAMBO.BasisSet_FunctionToCenter()) ...
    - centerNumFunctionsAO;
[indAOs, indAMBOs] = IndicesAOandAMBO(atomicNumbers, centerNumFunctionsAO, centerNumFunctionsAMBO);
overlapAOandAMBO = matpsi2AOandAMBO.Integrals_Overlap();
input.AOandAMBO = overlapAOandAMBO(indAOs, indAMBOs);

functionToCenterAOandAMBO = matpsi2AOandAMBO.BasisSet_FunctionToCenter();
input.functionToCenterAMBO = functionToCenterAOandAMBO(indAMBOs);

end

function centerNumFunctions = CenterNumFunctions(atomicNumbers, funcToCenter)
centerNumFunctions = zeros(length(atomicNumbers), 1);
for iatom = 1:length(atomicNumbers)
    centerNumFunctions(iatom) = sum(funcToCenter==iatom);
end
end

function [indAOs, indAMBOs] = IndicesAOandAMBO(atomicNumbers, centerNumFunctionsAO, centerNumFunctionsAMBO)
indAOs = [];
indAMBOs = [];
pointer = 1;
for iAtom = 1:length(atomicNumbers)
    numAOs = centerNumFunctionsAO(iAtom);
    numAMBOs = centerNumFunctionsAMBO(iAtom);
    indAOs = [indAOs pointer:pointer+numAOs-1]; %#ok
    pointer = pointer + numAOs;
    indAMBOs = [indAMBOs pointer:pointer+numAMBOs-1]; %#ok
    pointer = pointer + numAMBOs;
end
end
