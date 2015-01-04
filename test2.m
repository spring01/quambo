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
matpsi2 = MatPsi2(mol.MoleculeString,'6-31g*');

rhf = RHF(RHF.MatPsi2Interface(matpsi2));

rhf.SCF
rhf.DFSCF

tic
for i = 1:10
rhf.DFSCF;
end
toc

tic
for i = 1:10
rhf.SCF;
end
toc
