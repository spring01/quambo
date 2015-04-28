function [quambo, matpsi2AO] = MatPsi2Interface(molCart, basisSetInfo)

basisSetAO = basisSetInfo.basisSetAO;
basisSetAOandAMBO = [basisSetInfo.basisSetAO, '_and_', basisSetInfo.basisSetAMBO];
matpsi2AO = MatPsi2(molCart, basisSetAO);
matpsi2AO.SCF_RunRHF();
matpsi2AOandAMBO = MatPsi2(molCart, basisSetAOandAMBO, 0, 1, basisSetInfo.path);

input.numElectrons = matpsi2AO.Molecule_NumElectrons();
input.AOtoMO = matpsi2AO.SCF_OrbitalAlpha();

atomicNumbers = matpsi2AO.Molecule_AtomicNumbers();
centerNumFunctionsAO = ...
    CenterNumFunctions(atomicNumbers, matpsi2AO.BasisSet_FuncToCenter());
centerNumFunctionsAMBO = ...
    CenterNumFunctions(atomicNumbers, matpsi2AOandAMBO.BasisSet_FuncToCenter()) ...
    - centerNumFunctionsAO;
[indAOs, indAMBOs] = IndicesAOandAMBO(atomicNumbers, centerNumFunctionsAO, centerNumFunctionsAMBO);
overlapAOandAMBO = matpsi2AOandAMBO.Integrals_Overlap();
input.AOandAMBO = overlapAOandAMBO(indAOs, indAMBOs);

funcToCenterAOandAMBO = matpsi2AOandAMBO.BasisSet_FuncToCenter();
input.funcToCenterAMBO = funcToCenterAOandAMBO(indAMBOs);

quambo = QUAMBO(input);

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
