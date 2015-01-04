function hfEnergy = DFSCF(obj)
oeiVec = reshape(obj.kineticMat + sum(obj.potentialEachCoreMats,3), [], 1);
mn_Q = reshape(obj.mnATensor, [],size(obj.mnATensor,3)) * obj.invJHalfMetric;
m_nQ = reshape(mn_Q, size(obj.mnATensor,1),[]);

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
    
    temp = densKernelVec' * mn_Q;
    JVec = mn_Q * temp';
    
    temp = reshape( ...
        orbital(:,1:obj.numElectrons/2)' * m_nQ, ...
        obj.numElectrons/2, nbf,[]);
    temp = reshape(permute(temp, [2 1 3]), nbf,[]);
    KVec = reshape(temp * temp', [],1);
    
    fockVec = oeiVec + 2.*JVec - KVec; ... % - K
        
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
    disp('DFRHF.DoSCF(): Failed to converge.')
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

