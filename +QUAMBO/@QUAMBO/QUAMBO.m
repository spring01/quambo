classdef QUAMBO < handle
    
    properties (SetAccess = private)
        
        AOinQUAMBO;
        
        overlapQUAMBO;
        kineticQUAMBO;
        potentialEachCoreQUAMBO;
        twoElecIntegralsQUAMBO;
        
    end
    
    methods
        
        function obj = QUAMBO(properties)
            % Naming examples:
            % <AO|MO> : AOinMO
            % \sum_{MO}{<AO|MO><MO|AO>} : projMOinAO
            %
            % Diagonalize <AO|H|AO>, eigenkets are <AO|MO>
            %
            % MOinAO = AOinMO';
            % projMOinAO = AOinMO * AOinMO';
            
            % input section
            numOccMOs = properties.numElectrons / 2;
            numVirMOs = properties.numAOs - numOccMOs;
            numQUAMBOs = properties.numAMBOs;
            AOinMO = properties.AOinMO;
            AOinAMBO = properties.AOinAMBO;
            
            % diagonalize \sum_{AMBO}{<VirMO|AMBO><AMBO|VirMO>}
            % and take the significant virtual orbitals
            OccMOinAMBO = AOinMO(:, 1:numOccMOs)'*AOinAMBO; % <OccMO|AMBO> = \sum_{AO}{<OccMO|AO><AO|AMBO>}
            VirMOinAMBO = AOinMO(:, numOccMOs+1:end)'*AOinAMBO; % <VirMO|AMBO> = \sum_{AO}{<VirMO|AO><AO|AMBO>}
            projAMBOinVirMO = VirMOinAMBO * VirMOinAMBO'; % projector \sum{|AMBO><AMBO|} in |VirMO>
            [VirMOinWeightedVirMO, weights] = eig(projAMBOinVirMO);
            [~, order] = sort(diag(weights), 'descend');
            VirMOinSigVirMO = VirMOinWeightedVirMO(:, order(1:numQUAMBOs-numOccMOs));
            
            % maximize overlap
            projSigVirMOinVirMO = VirMOinSigVirMO * VirMOinSigVirMO';
            sumSqMOinAMBO = sum(OccMOinAMBO.^2) ... % \sum{OccMO}{<OccMO|AMBO>^2}
                + diag(VirMOinAMBO'*projSigVirMOinVirMO*VirMOinAMBO)'; % \sum{SigVirMO}{<AMBO|SigVirMO><SigVirMO|AMBO>}
            
            % calculate QUAMBOs
            OccMOinQUAMBO = repmat(sumSqMOinAMBO.^-0.5,numOccMOs,1) ...
                .* OccMOinAMBO; % <OccMO|QUAMBO>
            VirMOinQUAMBO = repmat(sumSqMOinAMBO.^-0.5,numVirMOs,1) ...
                .* (projSigVirMOinVirMO*VirMOinAMBO); % <VirMO|QUAMBO>
            MOinQUAMBO = [OccMOinQUAMBO; VirMOinQUAMBO]; % <MO|QUAMBO>
            
            % <QUAMBO|QUAMBO> and <AO|QUAMBO>
            obj.overlapQUAMBO = MOinQUAMBO' * MOinQUAMBO; % MOs are orthonormal
            obj.AOinQUAMBO = AOinMO * MOinQUAMBO;
            
            % transform other integrals
            obj.kineticQUAMBO = obj.AOinQUAMBO' * properties.kineticAO * obj.AOinQUAMBO;
            potentialEachCoreAO = properties.potentialEachCoreAO;
            obj.potentialEachCoreQUAMBO = zeros(numQUAMBOs,numQUAMBOs,size(potentialEachCoreAO, 3));
            for iAtom = 1:size(potentialEachCoreAO, 3)
                obj.potentialEachCoreQUAMBO(:,:,iAtom) = ...
                    obj.AOinQUAMBO' * potentialEachCoreAO(:,:,iAtom) * obj.AOinQUAMBO;
            end
            obj.twoElecIntegralsQUAMBO = QUAMBO.TransformTensor4(properties.twoElecIntegralsAO, obj.AOinQUAMBO);
        end
        
    end
    
    methods (Static)
        
        properties_ = UseMatPsi2(molStr, basisSetAO, basisSetAOandAMBO);
        
    end
    
    methods (Access = private)
        
        tensor4 = TransformTensor4(~, tensor4, trans);
        
    end
    
end