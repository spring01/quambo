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
        cartesian = mol.cartesian;
        theta1 = rand*200-100;
        rotMat1 = [cosd(theta1), -sind(theta1); sind(theta1), cosd(theta1)];
        cartesian(:, 2:4) = cartesian(:, 2:4) + rand;
        cartesian(:, 2:3) = (rotMat1*cartesian(:, 2:3)')';
        theta2 = rand*200-100;
        rotMat2 = [cosd(theta2), -sind(theta2); sind(theta2), cosd(theta2)];
        cartesian(:, 2:4) = cartesian(:, 2:4) + rand;
        cartesian(:, 3:4) = (rotMat2*cartesian(:, 3:4)')';
        theta3 = rand*200-100;
        rotMat3 = [cosd(theta3), -sind(theta3); sind(theta3), cosd(theta3)];
        cartesian(:, 2:4) = cartesian(:, 2:4) + rand;
        cartesian(:, [3,2]) = (rotMat3*cartesian(:, [3,2])')';
        cartesian(:, 2:4) = cartesian(:, 2:4) + rand;
        mol = Molecule(cartesian);
    end
    
    basisSetNames.basisSetAMBO = 'sto-3g-cartesian';
    basisSetNames.basisSetAO = '6-31gs';
    basisSetNames.auxiliaryBasisSet = 'cc-pvdz-jkfit';
    data1mol = RawDataFromOneMolecule(mol, basisSetNames);
    
    mp = MatPsi2(mol.MoleculeString, '6-31gs');
    
    comp{attempt} = data1mol.quambo.twoElecIntegralsQUAMBO;
%     comp{attempt} = data1mol.quambo.AOtoQUAMBO' * mp.Integrals_Overlap * data1mol.quambo.AOtoQUAMBO;
    
end

sum(sum(sum(sum(abs(abs(comp{1}) - abs(comp{2}))))))
sum(sum(sum(sum(abs(comp{1} - comp{2})))))

% sum(sum(sum(sum(abs(abs(comp{1}) - abs(comp{2}))))))
% sum(sum(sum(sum(abs(comp{1} - comp{2})))))

