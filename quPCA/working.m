load mols.mat

basisname = '6-31gs';

basisSetInfo.basisSetAO = '6-31gs';
basisSetInfo.basisSetAMBO = 'sto-3g-cartesian';
basisSetInfo.path = '/home/haichen/working/quambo/quPCA';

XTrain = zeros(0,25);
yTrain = zeros(0,15);
counter = 1;


for imol = 1:length(mols)
    molcart = mols{imol}.cartesian;
    
    matpsi = MatPsi2(molcart, basisname);
    matgdma = MatPsiGDMA(matpsi);
    matgdma.bigexp = inf;
    matpsi.SCF_RunRHF();
    psi4_occOrb = matpsi.SCF_OrbitalAlpha();
    psi4_occOrb = psi4_occOrb(:, 1:matpsi.Molecule_NumElectrons()/2);
    matgdma.RunGDMA(psi4_occOrb);
    
    matpsiAMBO = MatPsi2(molcart, basisSetInfo.basisSetAMBO, 0, 1, basisSetInfo.path);
    func2centerAMBO = matpsiAMBO.BasisSet_FuncToCenter;
    [quambo, matpsiAO] = QUAMBO.MatPsi2Interface(molcart, basisSetInfo);
    trans = quambo.AOtoQUAMBO;
    h1mat = trans' * matpsiAO.SCF_CoreHamiltonian() * trans;
    
    for iatom = 1:size(molcart, 1)
        if(molcart(iatom, 1) == 6) % is a carbon
            XTrain(counter, 1:25) = matgdma.multipoles(1:25, iatom);
            targetH1 = h1mat(func2centerAMBO == iatom, func2centerAMBO == iatom);
            yTrain(counter, 1:15) = targetH1(tril(true(5)));
            counter = counter + 1;
        end
    end
    
    disp(imol);
end



