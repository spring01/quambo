function properties = UseMatPsi2(molStr, basisSetAO, basisSetAOandAMBO)
matpsi2AO = MatPsi2(molStr, basisSetAO);

if(matpsi2AO.BasisSet_NumFunctions() > 150)
    throw(MException('QUAMBO:UseMatPsi2','Too many basis functions.'));
end

matpsi2AO.RHF_DoSCF();
matpsi2AOandAMBO = MatPsi2(molStr, basisSetAOandAMBO, 0, 1, [pwd(), '/+QUAMBO']);

properties.numElectrons = matpsi2AO.Molecule_NumElectrons();
properties.numAOs = matpsi2AO.BasisSet_NumFunctions();
properties.numAMBOs = matpsi2AOandAMBO.BasisSet_NumFunctions() - properties.numAOs;
properties.AOtoMO = matpsi2AO.RHF_Orbital();

atomicNumbers = matpsi2AO.Molecule_AtomicNumbers();
centerNumFunctionsAO = CenterNumFunctions(atomicNumbers, matpsi2AO.BasisSet_FunctionToCenter());
centerNumFunctionsAMBO = CenterNumFunctions(atomicNumbers, matpsi2AOandAMBO.BasisSet_FunctionToCenter()) ...
    - centerNumFunctionsAO;
[indAOs, indAMBOs] = RangeAOandAMBO(atomicNumbers, centerNumFunctionsAO, centerNumFunctionsAMBO);
overlapAOandAMBO = matpsi2AOandAMBO.Integrals_Overlap();
properties.AOandAMBO = overlapAOandAMBO(indAOs, indAMBOs);

properties.kineticAO = matpsi2AO.Integrals_Kinetic();
properties.potentialEachCoreAO = matpsi2AO.Integrals_PotentialEachCore();
properties.twoElecIntegralsAO = matpsi2AO.Integrals_AllTEIs();
end

function centerNumFunctions = CenterNumFunctions(atomicNumbers, funcToCenter)
centerNumFunctions = zeros(length(atomicNumbers), 1);
for iatom = 1:length(atomicNumbers)
    centerNumFunctions(iatom) = sum(funcToCenter==iatom);
end
end

function [indAOs, indAMBOs] = RangeAOandAMBO(atomicNumbers, centerNumFunctionsAO, centerNumFunctionsAMBO)
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
