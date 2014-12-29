import RawData.* FeatureLabel.*;

input = [...
    'O', char(10), ...
    'H    1    1.0', char(10), ...
    'H    1    1.0    2    120.0', char(10)];


mol = Molecule(input);

basisSetNames.minimalBasisSet = 'sto-3g-cartesian';
basisSetNames.largeBasisSet = '6-31gs';
data1mol = RawDataFromOneMolecule(mol, basisSetNames);
featLabel = FeatureLabelOverlap(data1mol, 1, 1);

