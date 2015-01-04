function hfEnergy = SCF(obj)
oeiVec = reshape(obj.kineticMat + sum(obj.potentialEachCoreMats,3), [], 1);
teiForCoulomb = reshape(obj.twoElecIntegrals, length(oeiVec), []);
teiForExchange = reshape(permute(obj.twoElecIntegrals, [1 3 2 4]), ...
    length(oeiVec), []);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densKernelVec = zeros(size(oeiVec));
elecEnergy = 0;

% diis adiis
cdiis = RawData.CDIIS(obj.overlapMat);
adiis = RawData.ADIIS(oeiVec);

fockVec = oeiVec;
fockSimVec = fockVec;
nbf = size(obj.kineticMat, 1);
for iter = 1:obj.maxSCFIter
    oldDensVec = densKernelVec;
    oldElecEnergy = elecEnergy;
    [densKernelVec, elecEnergy, orbital, orbitalEnergies] ...
        = DiagonalizeFock(reshape(fockSimVec, nbf, []), ...
        inv_S_Half, obj.numElectrons);
    elecEnergy = oeiVec'*densKernelVec + elecEnergy;
    
    if(sqrt(mean((densKernelVec - oldDensVec).^2)) < obj.RMSDensityThreshold ...
            && max(abs(densKernelVec - oldDensVec)) < obj.MaxDensityThreshold ...
            && abs(elecEnergy - oldElecEnergy) < obj.EnergyThreshold)
        break;
    end
    fockVec = oeiVec + ... % H
        2 .* (teiForCoulomb * densKernelVec) ... % + 2J
        - (teiForExchange * densKernelVec); ... % - K
        
    % diis extropolate Fock matrix
    cdiis.Push(fockVec, densKernelVec); % density must be idempotent
    adiis.Push(fockVec, densKernelVec); % Fock must be built from idempotent density
    if(cdiis.IAmBetter())
        fockSimVec = cdiis.Extrapolate();
    else
        fockSimVec = adiis.Interpolate();
    end
end

if(iter >= obj.maxSCFIter)
    disp('RHF.DoSCF(): Failed to converge.')
end

hfEnergy = elecEnergy + obj.nuclearRepulsionEnergy;

obj.finalFockMat = reshape(fockSimVec, nbf, []);
obj.finalDensKernelMat = reshape(densKernelVec, nbf, []);

obj.hfEnergy = hfEnergy;
obj.orbital = orbital;
obj.orbitalEnergies = orbitalEnergies;

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

