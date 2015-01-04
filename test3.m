import RawData.*;

input = [...
    'C', char(10), ...
    'C    1    1.54', char(10), ...
    'F    1    1.40    2    105.0', char(10), ...
    'H    1    1.08    2    105.0    3    120.0', char(10), ...
    'H    1    1.08    2    105.0    3   -120.0', char(10), ...
    'H    2    1.08    1    105.0    3     60.0', char(10), ...
    'H    2    1.08    1    105.0    3    -60.0', char(10), ...
    'H    2    1.08    1    105.0    3    180.0', char(10)];

mol = Molecule(input);
basisSetNames.basisSetAO = '6-31g*';
basisSetNames.basisSetAMBO = 'sto-3g-cartesian';
basisSetNames.auxiliaryBasisSet = 'cc-pvdz-jkfit';


matpsi2AO = MatPsi2(mol.MoleculeString,basisSetNames.basisSetAO);
matpsi2AO.JK_Initialize('DFJK', basisSetNames.auxiliaryBasisSet);
matpsi2AMBO = MatPsi2(mol.MoleculeString,basisSetNames.basisSetAMBO, 0,1,'./+RawData');
matpsi2AMBO.JK_Initialize('DFJK', basisSetNames.auxiliaryBasisSet);

quambo = QUAMBO(QUAMBO.UseMatPsi2(mol.MoleculeString, basisSetNames));

propQuambo = RHF.QUAMBOInterface(quambo);
propQuambo.invJHalfMetric = matpsi2AMBO.DFJK_InverseJHalfMetric();
propQuambo.nuclearRepulsionEnergy = matpsi2AMBO.Molecule_NuclearRepulsionEnergy();
propQuambo.numElectrons = matpsi2AMBO.Molecule_NumElectrons();


rhfquambo = RHF(propQuambo);
rhfao = RHF(RHF.MatPsi2Interface(matpsi2AO));

rhfquambo.SCF
rhfao.SCF
rhfquambo.DFSCF
rhfao.DFSCF

