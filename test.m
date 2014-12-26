import QUAMBO.*;

% molStr = [...
%     'H', char(10), ...
%     'F    1    1.0', char(10)];
molStr = [...
    'H  0.0  0.0  0.0', char(10), ...
    'F  0.0  0.0  1.0', char(10)];
molStr = [...
    'C', char(10), ...
    'C    1    1.54', char(10), ...
    'F    1    1.40    2    105.0', char(10), ...
    'H    1    1.07    2    105.0    3    120.0', char(10), ...
    'H    1    1.09    2    105.0    3   -120.0', char(10), ...
    'H    2    1.06    1    105.0    3     60.0', char(10), ...
    'H    2    1.08    1    105.0    3    -60.0', char(10), ...
    'H    2    1.04    1    105.0    3    180.0', char(10)];

basisSet = '6-31gs';
mpsi2 = MatPsi2(molStr, basisSet);


basisSetAO = '6-31gs';
basisSetAOandAMBO = '6-31gs_and_sto-3g-cartesian';

quambo = QUAMBO(QUAMBO.UseMatPsi2(molStr, basisSetAO, basisSetAOandAMBO));


properties.overlapMat = quambo.overlapQUAMBO;
properties.kineticMat = quambo.kineticQUAMBO;
properties.corePotentialMat = sum(quambo.potentialEachCoreQUAMBO, 3);
properties.twoElecIntegrals = quambo.twoElecIntegralsQUAMBO;
properties.numElectrons = mpsi2.Molecule_NumElectrons();
properties.nuclearRepulsionEnergy = mpsi2.Molecule_NuclearRepulsionEnergy();

rhf = RHF(properties);
[orbital, orbitalEnergies, hfEnergy, iter] = rhf.SCF();
hfEnergy
mpsi2.RHF_DoSCF()

