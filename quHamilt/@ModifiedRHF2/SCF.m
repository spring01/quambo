function [hfEnergy, orbital, orbitalEnergies, iter] = SCF(obj)
oeiVec = reshape(obj.kineticMat .* obj.kineticMatModifier ...
    + sum(obj.corePotentialMats .* obj.corePotentialMatsModifier,3), [], 1);
teiForCoulomb = reshape(obj.twoElecIntegrals .* obj.teiModifier, length(oeiVec), []);
teiForExchange = reshape(permute(obj.twoElecIntegrals .* obj.teiModifier, [1 3 2 4]), ...
    length(oeiVec), []);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat .* obj.overlapMatModifier);

densVec = zeros(size(oeiVec));
elecEnergy = 0;

% diis adiis
cdiis = CDIIS(obj.overlapMat);
adiis = ADIIS(oeiVec);

fockVec = oeiVec;
fockSimVec = fockVec;
for iter = 1:obj.maxSCFIter
    oldDensVec = densVec;
    oldElecEnergy = elecEnergy;
    [densVec, elecEnergy, orbital, orbitalEnergies] ...
        = DiagonalizeFock(reshape(fockSimVec, sqrt(length(fockSimVec)), []), ...
        inv_S_Half, obj.numElectrons);
    elecEnergy = oeiVec'*densVec + elecEnergy;
    
    if(sqrt(mean((densVec - oldDensVec).^2)) < obj.RMSDensityThreshold ...
            && max(abs(densVec - oldDensVec)) < obj.MaxDensityThreshold ...
            && abs(elecEnergy - oldElecEnergy) < obj.EnergyThreshold)
        break;
    end
    fockVec = oeiVec + ... H
        2 .* (teiForCoulomb * densVec) ... % +2J
        - (teiForExchange * densVec); ... % -K
        
    % diis extropolate Fock matrix
    cdiis.Push(fockVec, densVec); % density must be idempotent
    adiis.Push(fockVec, densVec); % Fock must be built from idempotent density
    if(cdiis.IAmBetter())
        fockSimVec = cdiis.Extrapolate();
    else
        fockSimVec = adiis.Interpolate();
    end
end
hfEnergy = elecEnergy + obj.nucRepEnergy;

obj.hfEnergy = hfEnergy;
obj.orbital = orbital;
obj.orbitalEnergies = orbitalEnergies;
obj.coulombMat = reshape(teiForCoulomb * densVec, size(obj.overlapMat));
obj.exchangeMat = reshape(teiForExchange * densVec, size(obj.overlapMat));

end


function [densityVec, elecEnergy, orbital, orbitalEnergies] = DiagonalizeFock(fockMat, inv_S_Half, numElectrons)
[orbitalOtho, orbitalEnergies] = eig(inv_S_Half*fockMat*inv_S_Half);
[orbitalEnergies, ascend_order] = sort(diag(orbitalEnergies));
orbital = inv_S_Half * orbitalOtho(:, ascend_order);
densityVec = reshape( ...
    orbital(:, 1:numElectrons/2) * orbital(:, 1:numElectrons/2)', ...
    [], 1);
elecEnergy = sum(orbitalEnergies(1:numElectrons/2));
end

