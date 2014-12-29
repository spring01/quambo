basisSet = '6-31gs';
mpsi2 = MatPsi2(mol.MoleculeString, basisSet);


basisSetAO = '6-31gs';
basisSetAOandAMBO = '6-31gs_and_sto-3g-cartesian';

quambo = QUAMBO(QUAMBO.UseMatPsi2(mol.MoleculeString, basisSetAO, basisSetAOandAMBO, '/home/haichen/working/quambo/+RawData'));


mpsi2MBS = MatPsi2(mol.MoleculeString, 'sto-3g-cartesian', 0, 1, '/home/haichen/working/quambo/+RawData');

trans = quambo.AMBOtoQUAMBO;
properties.overlapMat = trans' * mpsi2MBS.Integrals_Overlap * trans;
properties.kineticMat = trans' * mpsi2MBS.Integrals_Kinetic * trans;
properties.potentialMat = trans' * mpsi2MBS.Integrals_Potential * trans;
properties.twoElecIntegrals = TransformTensor4(mpsi2MBS.Integrals_AllTEIs, trans);
properties.numElectrons = mpsi2MBS.Molecule_NumElectrons();
properties.nuclearRepulsionEnergy = mpsi2MBS.Molecule_NuclearRepulsionEnergy();

rhf = RHF(properties);
hfEnergy = rhf.DoSCF();
hfEnergy
mpsi2.RHF_DoSCF()