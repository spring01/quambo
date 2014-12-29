import RawData.*;

input = [...
    'O', char(10), ...
    'H    1    1.0', char(10), ...
    'H    1    1.0    2    120.0', char(10)];
input = [...
    'H  0.0  0.0  0.0', char(10), ...
    'F  0.0  0.0  1.0', char(10)];
input = [...
    'C', char(10), ...
    'C    1    1.54', char(10), ...
    'F    1    1.40    2    105.0', char(10), ...
    'H    1    1.07    2    105.0    3    120.0', char(10), ...
    'H    1    1.09    2    105.0    3   -120.0', char(10), ...
    'H    2    1.06    1    105.0    3     60.0', char(10), ...
    'H    2    1.08    1    105.0    3    -60.0', char(10), ...
    'H    2    1.04    1    105.0    3    180.0', char(10)];

input = [...
    'C', ',', ...
    'C    1    1.54', char(10), ...
    'F    1    1.40    2    105.0', ';', ...
    'H    1    1.07    2    105.0    3    120.0', char(10), ...
    'H    1    1.09    2    105.0    3   -120.0', char(10), ...
    'H    2    1.06    1    105.0    3     60.0', ':', ...
    'H    2    1.08    1    105.0    3    -60.0', '}', ...
    'H    2    1.04    1    105.0    3    180.0', char(10)];

% input = [ ...
%     6 0 0.0 0 0.0 0  0.0; ...
%     1 1 1.0 0 0.0 0  0.0; ...
%     1 1 1.0 2 105 0  0.0; ...
%     1 1 1.0 2 105 3  120; ...
%     1 1 1.0 2 105 3 -120;];
% 
% input = [...
%     'C', char(10), ...
%     'C    1    1.4', char(10), ...
%     'C    2    1.4    1    120.0', char(10), ...
%     'C    3    1.4    2    120.0    1      0.0', char(10), ...
%     'C    4    1.4    3    120.0    2      0.0', char(10), ...
%     'C    5    1.4    4    120.0    3      0.0', char(10), ...
%     'H    1    1.0    6    120.0    5    180.0', char(10), ...
%     'H    2    1.0    1    120.0    6    180.0', char(10), ...
%     'H    3    1.0    2    120.0    1    180.0', char(10), ...
%     'H    4    1.0    3    120.0    2    180.0', char(10), ...
%     'H    5    1.0    4    120.0    3    180.0', char(10), ...
%     'H    6    1.0    5    120.0    4    180.0'];

mol = Molecule(input);

basisSets.minimalBasisSet = 'sto-3g-cartesian';
basisSets.largeBasisSet = '6-31gs';
data1mol = RawDataFromOneMolecule(mol, basisSets, '/home/haichen/working/quambo/+RawData');

% basisSet = '6-31gs';
% mpsi2 = MatPsi2(mol.MoleculeString, basisSet);
% 
% 
% basisSetAO = '6-31gs';
% basisSetAOandAMBO = '6-31gs_and_sto-3g-cartesian';
% 
% quambo = QUAMBO(QUAMBO.UseMatPsi2(mol.MoleculeString, basisSetAO, basisSetAOandAMBO));
% 
% 
% properties.overlapMat = quambo.overlapQUAMBO;
% properties.kineticMat = quambo.kineticQUAMBO;
% properties.corePotentialMat = sum(quambo.potentialEachCoreQUAMBO, 3);
% properties.twoElecIntegrals = quambo.twoElecIntegralsQUAMBO;
% properties.numElectrons = mpsi2.Molecule_NumElectrons();
% properties.nuclearRepulsionEnergy = mpsi2.Molecule_NuclearRepulsionEnergy();
% 
% rhf = RHF(properties);
% [orbital, orbitalEnergies, hfEnergy, iter] = rhf.SCF();
% hfEnergy
% mpsi2.RHF_DoSCF()
