classdef QUAMBO < handle
    
    properties (SetAccess = private)
        
        AOtoQUAMBO;
        
        overlapQUAMBO;
        kineticQUAMBO;
        potentialEachCoreQUAMBO;
        twoElecIntegralsQUAMBO;
        
    end
    
    methods
        
        function obj = QUAMBO(properties)
            
            funcToCenter = properties.functionToCenterAMBO;
            hamiltonianAMBO = properties.hamiltonianAMBO;
            AOtoMO = properties.AOtoMO;
            AOandAMBO = properties.AOandAMBO;
            numElectrons = properties.numElectrons;
            
            overlapAO = properties.overlapAO;
            kineticAO = properties.kineticAO;
            potentialEachCoreAO = properties.potentialEachCoreAO;
            
            MOtoQUAMBO = obj.FindMOtoQUAMBO(AOtoMO, AOandAMBO, numElectrons);
            MOtoQUAMBO = obj.SymmOrthoNShells(MOtoQUAMBO, funcToCenter);
            AOtoQUAMBO_ = AOtoMO * MOtoQUAMBO;                               % <AO\MO><MO\QUAMBO> = <AO\QUAMBO>
            
            tempHamiltQUAMBO ...
                = AOtoQUAMBO_' ...
                * (kineticAO+sum(potentialEachCoreAO,3)) ...
                * AOtoQUAMBO_;
            invarTransQUAMBO = obj.InvariantTransform( ...
                tempHamiltQUAMBO, hamiltonianAMBO, funcToCenter);
            invarTransAMBO = obj.InvariantTransform( ...
                hamiltonianAMBO, hamiltonianAMBO, funcToCenter);
            
            AOtoQUAMBO_ = AOtoQUAMBO_ * invarTransQUAMBO * invarTransAMBO';
            
            % integrals
            obj.overlapQUAMBO = AOtoQUAMBO_' * overlapAO * AOtoQUAMBO_;
            obj.kineticQUAMBO = AOtoQUAMBO_' * kineticAO * AOtoQUAMBO_;
            obj.potentialEachCoreQUAMBO = zeros( ...
                [size(obj.kineticQUAMBO), ...
                size(potentialEachCoreAO, 3)]);
            for iAtom = 1:size(potentialEachCoreAO, 3)
                obj.potentialEachCoreQUAMBO(:,:,iAtom) ...
                    = AOtoQUAMBO_' ...
                    * potentialEachCoreAO(:,:,iAtom) ...
                    * AOtoQUAMBO_;
            end
            obj.twoElecIntegralsQUAMBO = obj.TransformTensor4( ...
                properties.twoElecIntegralsAO, AOtoQUAMBO_);
            
            obj.AOtoQUAMBO = AOtoQUAMBO_;
        end
        
    end
    
    methods (Access = private)
        
        function MOtoQUAMBO = FindMOtoQUAMBO(~, AOtoMO, AOandAMBO, numElectrons)
            numQUAMBOs = size(AOandAMBO, 2);
            numOccMOs = numElectrons / 2;
            numVirMOs = size(AOandAMBO, 1) - numOccMOs;
            
            % diagonalize \sum{|AMBO><AMBO|} in |VirMO> and take the significant virtual orbitals
            OccMOtoAMBO = AOtoMO(:, 1:numOccMOs)'*AOandAMBO;                % <OccMO/AO><AO|AMBO> = <OccMO|AMBO> = <OccMO\AMBO>
            VirMOtoAMBO = AOtoMO(:, numOccMOs+1:end)'*AOandAMBO;            % <VirMO/AO><AO|AMBO> = <VirMO|AMBO> = <VirMO\AMBO>
            [VirMOtoWeightedVirMO, weights] = ...
                eig(VirMOtoAMBO*VirMOtoAMBO');                              % <VirMO|AMBO><AMBO|VirMO>
            [~, order] = sort(diag(weights), 'descend');
            VirMOtoSigVirMO = ...
                VirMOtoWeightedVirMO(:, order(1:numQUAMBOs-numOccMOs));     % <VirMO\SigVirMO>
            
            % maximize overlap
            SigVirMOtoAMBO = VirMOtoSigVirMO'*VirMOtoAMBO;                  % <SigVirMO/VirMO><VirMO|AMBO> = <SigVirMO|AMBO> = <SigVirMO\AMBO>
            sumMOtoAMBOsq = sum([OccMOtoAMBO;SigVirMOtoAMBO].^2);
            
            % QUAMBOs
            OccMOtoQUAMBO = repmat(sumMOtoAMBOsq.^-0.5,numOccMOs,1) ...     % <OccMO\QUAMBO> ~ <OccMO\AMBO>
                .* OccMOtoAMBO;
            VirMOtoQUAMBO = repmat(sumMOtoAMBOsq.^-0.5,numVirMOs,1) ...
                .* (VirMOtoSigVirMO*SigVirMOtoAMBO);                        % <VirMO\QUAMBO> ~ <VirMO\SigVirMO><SigVirMO|AMBO>
            MOtoQUAMBO = [OccMOtoQUAMBO; VirMOtoQUAMBO];
        end
        
        function MOtoQUAMBO = SymmOrthoNShells(~, MOtoQUAMBO, funcToCenter)
            tempOverlap = MOtoQUAMBO' * MOtoQUAMBO;
            for iAtom = unique(funcToCenter)
                funcs = find(funcToCenter==iAtom);
                if(length(funcs) > 1) % orthogonalize 2s+2p
                    MOtoQUAMBO(:,funcs(2:5)) = MOtoQUAMBO(:,funcs(2:5)) ...
                        / sqrtm(tempOverlap(funcs(2:5),funcs(2:5)));
                end
                if(length(funcs) > 5) % orthogonalize 3s+3p
                    MOtoQUAMBO(:,funcs(6:9)) = MOtoQUAMBO(:,funcs(6:9)) ...
                        / sqrtm(tempOverlap(funcs(6:9),funcs(6:9)));
                end
                if(length(funcs) > 9)
                    throw(MException('QUAMBO:SymmOrthoNShells' , ...
                        'Cannot handle 4th(+) row elements.'));
                end
            end
        end
        
        function trans = ...
                InvariantTransform(~, hamilt, hamiltForVecs, funcToCenter)
            trans = eye(size(hamiltForVecs));
            for iAtom = unique(funcToCenter)
                funcs = find(funcToCenter==iAtom);
                if(length(funcs) > 1) % orthogonalize 2s+2p
                    [vec, ~] = eig(hamiltForVecs(funcs(3:5),funcs(3:5)));
                    trans(:,funcs(3:5)) = trans(:,funcs(3:5)) * vec;
                end
                if(length(funcs) > 5) % orthogonalize 3s+3p
                    [vec, ~] = eig(hamiltForVecs(funcs(7:9),funcs(7:9)));
                    trans(:,funcs(7:9)) = trans(:,funcs(7:9)) * vec;
                end
                if(length(funcs) > 9)
                    throw(MException('QUAMBO:SymmOrthoNShells' , ...
                        'Cannot handle 4th(+) row elements.'));
                end
            end
            invarHamilt = trans' * hamilt * trans;
            
            firstHeavyZFound = false;
            firstHydrogenFound = false;
            for iAtom = unique(funcToCenter)
                funcs = find(funcToCenter==iAtom);
                if(length(funcs) > 1) % orthogonalize 2s+2p
                    if(~firstHeavyZFound)
                        firstHeavyZFound = true;
                        firstHeavyZOrbital = funcs(5);
                        firstHeavyAtomIndex = iAtom;
                    end
                    if(invarHamilt(funcs(3), funcs(2)) > 0)
                        trans(:,funcs(3)) = -trans(:,funcs(3));
                    end
                    if(invarHamilt(funcs(4), funcs(2)) > 0)
                        trans(:,funcs(4)) = -trans(:,funcs(4));
                    end
                else
                    if(~firstHydrogenFound)
                        firstHydrogenFound = true;
                        firstHydrogenOrbital = funcs;
                    end
                end
                if(length(funcs) > 5) % orthogonalize 3s+3p
                    if(invarHamilt(funcs(7), funcs(6)) > 0)
                        trans(:,funcs(7)) = -trans(:,funcs(7));
                    end
                    if(invarHamilt(funcs(8), funcs(6)) > 0)
                        trans(:,funcs(8)) = -trans(:,funcs(8));
                    end
                end
            end
            
            if(invarHamilt(firstHeavyZOrbital, firstHydrogenOrbital) > 0)
                trans(:,firstHeavyZOrbital) = -trans(:,firstHeavyZOrbital);
            end
            invarHamilt = trans' * hamilt * trans;
            
            refHeavyZOrbital = firstHeavyZOrbital;
            allAtoms = unique(funcToCenter);
            for iAtom = allAtoms(allAtoms~=firstHeavyAtomIndex)
                funcs = find(funcToCenter==iAtom);
                if(length(funcs) > 1) % orthogonalize 2s+2p
                    if(invarHamilt(funcs(5), refHeavyZOrbital) > 0)
                        trans(:,funcs(5)) = -trans(:,funcs(5));
                        invarHamilt = trans' * hamilt * trans;
                    end
                    refHeavyZOrbital = funcs(5);
                end
                if(length(funcs) > 5) % orthogonalize 3s+3p
                    if(invarHamilt(funcs(9), refHeavyZOrbital) > 0)
                        trans(:,funcs(9)) = -trans(:,funcs(9));
                        invarHamilt = trans' * hamilt * trans;
                    end
                    refHeavyZOrbital = funcs(9);
                end
            end
        end
        
        function tensor4 = TransformTensor4(~, tensor4, trans)
            nbf1 = size(trans, 1);
            nbf2 = size(trans, 2);
            
            % pqrs -> iqrs
            tensor4 = reshape( ...
                trans' * reshape(tensor4, nbf1,[]) ...
                , nbf2,nbf1,nbf1,nbf1);
            
            % iqrs -> ijrs
            tensor4 = permute(tensor4,[2 1 3 4]);
            tensor4 = trans' * reshape(tensor4, nbf1,[]);
            tensor4 = permute( ...
                reshape(tensor4, nbf2,nbf2,nbf1,nbf1) ...
                , [2 1 3 4]);
            
            % ijrs -> ijks
            tensor4 = permute(tensor4,[1 2 4 3]);
            tensor4 = reshape(tensor4, [],nbf1) * trans;
            tensor4 = permute( ...
                reshape(tensor4, nbf2,nbf2,nbf1,nbf2) ...
                , [1 2 4 3]);
            
            % ijks -> ijkl
            tensor4 = reshape(tensor4, [],nbf1) * trans;
            tensor4 = reshape(tensor4, nbf2,nbf2,nbf2,nbf2);
        end

    end
    
    methods (Static)
        
        properties_ = UseMatPsi2(molStr, basisSetAO, basisSetAOandAMBO);
        
    end
    
end
