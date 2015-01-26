import RawData.*;
import FeatureLabel.*;


inputMol = [...
    'C', char(10), ...
    'C    1    1.54', char(10), ...
    'H    1    1.08    2    105.0', char(10), ...
    'H    1    1.08    2    105.0    3    120.0', char(10), ...
    'H    1    1.08    2    105.0    3   -120.0', char(10), ...
    'H    2    1.08    1    105.0    3     60.0', char(10), ...
    'H    2    1.08    1    105.0    3    -60.0', char(10), ...
    'F    2    1.48    1    105.0    3    180.0', char(10)];

mol = Molecule(inputMol);

basisSetNames.basisSetAO = '6-31gs';
basisSetNames.basisSetAMBO = 'sto-3g-cartesian';

rawData1mol = RawData1Mol(mol, basisSetNames);



mlh = rawData1mol.mlHueckelTarget;
matpsiAO = MatPsi2(mol.MoleculeString(), basisSetNames.basisSetAO);
matpsiAO.RHF_DoSCF();

flo1 = FeatureLabelOverlapOffDiag(rawData1mol, 1,2);

flhd1 = FeatureLabelCoreHamiltDiag(rawData1mol, 1);
flho1 = FeatureLabelCoreHamiltOffDiag(rawData1mol, 1,2);
flfd1 = FeatureLabelFockDiag(rawData1mol, 1);
flfo1 = FeatureLabelFockOffDiag(rawData1mol, 1,2);

