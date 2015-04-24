
inputMol = [...
    'C', char(10), ...
    'C    1    1.54', char(10), ...
    'H    1    1.08    2    105.0', char(10), ...
    'H    1    1.08    2    105.0    3    120.0', char(10), ...
    'H    1    1.08    2    105.0    3   -120.0', char(10), ...
    'H    2    1.08    1    105.0    3     60.0', char(10), ...
    'H    2    1.08    1    105.0    3    -60.0', char(10), ...
    'H    2    1.08    1    105.0    3    180.0', char(10)];

mol = Molecule(inputMol);

basisSetInfo.basisSetAO = '6-31gs';
basisSetInfo.basisSetAMBO = 'sto-3g-cartesian';
basisSetInfo.path = './';

mphl = MatPsi2(mol.cartesian, basisSetInfo.basisSetAO);
rhfhl = RHF.MatPsi2Interface(mphl);



quambo = QUAMBO.MatPsi2Interface(mol, basisSetInfo);
rhfqu = quambo.NewRHF(rhfhl);
rhfqu.SCF();

mpll = MatPsi2(mol.cartesian, basisSetInfo.basisSetAMBO, 0, 1, basisSetInfo.path);

modrhf = ModifiedRHF(mpll, rhfqu);
modrhf.SCF
