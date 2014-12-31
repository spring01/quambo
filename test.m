import RawData.* FeatureLabel.*;

input = [...
    'O    0              0         0', char(10), ...
    'H    0              0    1.0000', char(10), ...
    'H    0.8660         0   -0.5000', char(10)];

input = [...
    'C', char(10), ...
    'C    1    1.54', char(10), ...
    'F    1    1.40    2    105.0', char(10), ...
    'H    1    1.08    2    105.0    3    120.0', char(10), ...
    'H    1    1.08    2    105.0    3   -120.0', char(10), ...
    'H    2    1.08    1    105.0    3     60.0', char(10), ...
    'H    2    1.08    1    105.0    3    -60.0', char(10), ...
    'H    2    1.08    1    105.0    3    180.0', char(10)];

% input = [...
%     'C', char(10), ...
%     'C    1    1.44', char(10), ...
%     'F    1    1.40    2    120.0', char(10), ...
%     'H    1    1.08    2    120.0    3    180.0', char(10), ...
%     'H    2    1.08    1    120.0    3      0.0', char(10), ...
%     'H    2    1.08    1    120.0    3    180.0', char(10)];

% input = [...
%     'C', char(10), ...
%     'C    1    1.24', char(10), ...
%     'F    1    1.40    2    179.0', char(10), ...
%     'H    2    1.07    1    180.0    3    180.0'];

for attempt = 1:2
    
    mol = Molecule(input);
    
    if(attempt == 2)
        theta1 = rand*200-100;
        rotMat1 = [cosd(theta1), -sind(theta1); sind(theta1), cosd(theta1)];
        mol.cartesian(:, 2:4) = mol.cartesian(:, 2:4) + rand;
        mol.cartesian(:, 2:3) = (rotMat1*mol.cartesian(:, 2:3)')';
        theta2 = rand*200-100;
        rotMat2 = [cosd(theta2), -sind(theta2); sind(theta2), cosd(theta2)];
        mol.cartesian(:, 2:4) = mol.cartesian(:, 2:4) + rand;
        mol.cartesian(:, 3:4) = (rotMat2*mol.cartesian(:, 3:4)')';
        theta3 = rand*200-100;
        rotMat3 = [cosd(theta3), -sind(theta3); sind(theta3), cosd(theta3)];
        mol.cartesian(:, 2:4) = mol.cartesian(:, 2:4) + rand;
        mol.cartesian(:, [3,2]) = (rotMat3*mol.cartesian(:, [3,2])')';
        mol.cartesian(:, 2:4) = mol.cartesian(:, 2:4) + rand;
    end
    
    basisSetNames.minimalBasisSet = 'sto-3g-cartesian';
    basisSetNames.largeBasisSet = '6-31gs';
    data1mol = RawDataFromOneMolecule(mol, basisSetNames);
    
    mp = MatPsi2(mol.MoleculeString, '6-31gs');
    
    AOtoQUAMBO = data1mol.quambo.AOtoQUAMBO;
    
    inds1 = 2:5;
    inds2 = 7:10;
    inds3 = 12:15;
    invshalf1 = eye(length(inds1))/sqrtm(data1mol.quambo.overlapQUAMBO(inds1,inds1));
    invshalf2 = eye(length(inds2))/sqrtm(data1mol.quambo.overlapQUAMBO(inds2,inds2));
    invshalf3 = eye(length(inds3))/sqrtm(data1mol.quambo.overlapQUAMBO(inds3,inds3));
    AOtoQUAMBO(:,inds1) = AOtoQUAMBO(:,inds1)*invshalf1;
    AOtoQUAMBO(:,inds2) = AOtoQUAMBO(:,inds2)*invshalf2;
    AOtoQUAMBO(:,inds3) = AOtoQUAMBO(:,inds3)*invshalf3;
    
    
    tempover = AOtoQUAMBO'*(mp.Integrals_Overlap)*AOtoQUAMBO;
    if(tempover(2,1) < 0)
        AOtoQUAMBO(:,2) = -AOtoQUAMBO(:,2);
    end
    tempover = AOtoQUAMBO'*(mp.Integrals_Overlap)*AOtoQUAMBO;
    for i = [6 7 11 12 16 17 18 19 20]
        if(tempover(i,2) < 0)
            AOtoQUAMBO(:,i) = -AOtoQUAMBO(:,i);
        end
    end
    
    inds1 = 3:5;
    inds2 = 8:10;
    inds3 = 13:15;
    mpmbs = MatPsi2(mol.MoleculeString, 'sto-3g-cartesian',0,1,[pwd, '/+RawData']);
    tempkin = (mpmbs.Integrals_Kinetic + mpmbs.Integrals_Potential);
    [vec1, val1] = eig(tempkin(inds1,inds1));
    [vec2, val2] = eig(tempkin(inds2,inds2));
    [vec3, val3] = eig(tempkin(inds3,inds3));
    AOtoQUAMBO(:,inds1) = AOtoQUAMBO(:,inds1)*vec1;
    AOtoQUAMBO(:,inds2) = AOtoQUAMBO(:,inds2)*vec2;
    AOtoQUAMBO(:,inds3) = AOtoQUAMBO(:,inds3)*vec3;
    
    tempover = AOtoQUAMBO'*(mp.Integrals_Kinetic + mp.Integrals_Potential)*AOtoQUAMBO;
    for i = 1:length(inds1(1:2))
        index = inds1(i);
        if(tempover(index,inds1(1)-1) > 0)
            AOtoQUAMBO(:,index) = -AOtoQUAMBO(:,index);
            vec1(:,i) = -vec1(:,i);
        end
    end
    for i = 1:length(inds2(1:2))
        index = inds2(i);
        if(tempover(index,inds2(1)-1) > 0)
            AOtoQUAMBO(:,index) = -AOtoQUAMBO(:,index);
            vec2(:,i) = -vec2(:,i);
        end
    end
    for i = 1:length(inds3(1:2))
        index = inds3(i);
        if(tempover(index,inds3(1)-1) > 0)
            AOtoQUAMBO(:,index) = -AOtoQUAMBO(:,index);
            vec3(:,i) = -vec3(:,i);
        end
    end
    % C1,2p3
    if(tempover(5,16) > 0)
        AOtoQUAMBO(:,5) = -AOtoQUAMBO(:,5);
        vec1(:,3) = -vec1(:,3);
        tempover = AOtoQUAMBO'*(mp.Integrals_Kinetic + mp.Integrals_Potential)*AOtoQUAMBO;
    end
    % C2,2p3
    if(tempover(10,5) > 0)
        AOtoQUAMBO(:,10) = -AOtoQUAMBO(:,10);
        vec2(:,3) = -vec2(:,3);
    end
    % F1,2p3
    if(tempover(15,5) > 0)
        AOtoQUAMBO(:,15) = -AOtoQUAMBO(:,15);
        vec3(:,3) = -vec3(:,3);
    end
    
%     % now we have "invariant" integrals; let's transform them back
%     AOtoQUAMBO(:,inds1) = AOtoQUAMBO(:,inds1)*vec1';
%     AOtoQUAMBO(:,inds2) = AOtoQUAMBO(:,inds2)*vec2';
%     AOtoQUAMBO(:,inds3) = AOtoQUAMBO(:,inds3)*vec3';
    
    
    
    properties.kineticMat = AOtoQUAMBO'*mp.Integrals_Kinetic*AOtoQUAMBO;
    temp = mp.Integrals_PotentialEachCore;
    properties.potentialEachCoreMats = zeros([size(properties.kineticMat), mp.Molecule_NumAtoms]);
    for i = 1:mp.Molecule_NumAtoms
        properties.potentialEachCoreMats(:,:,i) = AOtoQUAMBO'*temp(:,:,i)*AOtoQUAMBO;
    end
    properties.twoElecIntegrals = TransformTensor4(mp.Integrals_AllTEIs, AOtoQUAMBO);
    properties.overlapMat = AOtoQUAMBO'*mp.Integrals_Overlap*AOtoQUAMBO;
    properties.numElectrons = mp.Molecule_NumElectrons;
    properties.nuclearRepulsionEnergy = mp.Molecule_NuclearRepulsionEnergy;
    
    precision = 1e8;
    properties.kineticMat = round(properties.kineticMat.*precision)./precision;
    properties.potentialEachCoreMats = round(properties.potentialEachCoreMats.*precision)./precision;
    properties.twoElecIntegrals = round(properties.twoElecIntegrals.*precision)./precision;
    properties.overlapMat = round(properties.overlapMat.*precision)./precision;
    
    rhftest = RHF(properties);
    627*(rhftest.DoSCF - mp.RHF_DoSCF)
    
    comp{attempt} = properties.overlapMat;
%     comp{attempt} = properties.kineticMat;
    comp{attempt} = properties.potentialEachCoreMats;
    comp{attempt} = properties.twoElecIntegrals;
    
end

sum(sum(sum(sum(abs(abs(comp{1}) - abs(comp{2}))))))
sum(sum(sum(sum(abs(comp{1} - comp{2})))))

