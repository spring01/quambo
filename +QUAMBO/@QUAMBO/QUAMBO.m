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
            % input section
            numOccMOs = properties.numElectrons / 2;
            numVirMOs = properties.numAOs - numOccMOs;
            numQUAMBOs = properties.numAMBOs;
            AOtoMO = properties.AOtoMO;
            AOandAMBO = properties.AOandAMBO;
            
            % diagonalize \sum{|AMBO><AMBO|} in |VirMO> and take the significant virtual orbitals
            OccMOtoAMBO = AOtoMO(:, 1:numOccMOs)'*AOandAMBO; % <OccMO/AO><AO|AMBO> = <OccMO|AMBO> = <OccMO\AMBO>
            VirMOtoAMBO = AOtoMO(:, numOccMOs+1:end)'*AOandAMBO; % <VirMO/AO><AO|AMBO> = <VirMO|AMBO> = <VirMO\AMBO>
            [VirMOtoWeightedVirMO, weights] = eig(VirMOtoAMBO*VirMOtoAMBO'); % <VirMO|AMBO><AMBO|VirMO>
            [~, order] = sort(diag(weights), 'descend');
            VirMOtoSigVirMO = VirMOtoWeightedVirMO(:, order(1:numQUAMBOs-numOccMOs)); % <VirMO\SigVirMO>
            
            % maximize overlap
            SigVirMOtoAMBO = VirMOtoSigVirMO'*VirMOtoAMBO; % <SigVirMO/VirMO><VirMO|AMBO> = <SigVirMO|AMBO> = <SigVirMO\AMBO>
            sumMOtoAMBOsq = sum([OccMOtoAMBO;SigVirMOtoAMBO].^2);
            
            % QUAMBOs
            OccMOtoQUAMBO = repmat(sumMOtoAMBOsq.^-0.5,numOccMOs,1) .* OccMOtoAMBO; % <OccMO\QUAMBO> ~ <OccMO\AMBO>
            VirMOtoQUAMBO = repmat(sumMOtoAMBOsq.^-0.5,numVirMOs,1) .* (VirMOtoSigVirMO*SigVirMOtoAMBO); % <VirMO\QUAMBO> ~ <VirMO\SigVirMO><SigVirMO|AMBO>
            MOtoQUAMBO = [OccMOtoQUAMBO; VirMOtoQUAMBO];
            obj.AOtoQUAMBO = AOtoMO * MOtoQUAMBO; % <AO\MO><MO\QUAMBO> = <AO\QUAMBO>
            
            % integrals
            obj.overlapQUAMBO = MOtoQUAMBO' * MOtoQUAMBO;
            obj.kineticQUAMBO = obj.AOtoQUAMBO' * properties.kineticAO * obj.AOtoQUAMBO;
            obj.potentialEachCoreQUAMBO = zeros(numQUAMBOs,numQUAMBOs,size(properties.potentialEachCoreAO, 3));
            for iAtom = 1:size(properties.potentialEachCoreAO, 3)
                obj.potentialEachCoreQUAMBO(:,:,iAtom) = ...
                    obj.AOtoQUAMBO' * properties.potentialEachCoreAO(:,:,iAtom) * obj.AOtoQUAMBO;
            end
            obj.twoElecIntegralsQUAMBO = ...
                QUAMBO.TransformTensor4(properties.twoElecIntegralsAO, obj.AOtoQUAMBO);
        end
        
    end
    
    methods (Static)
        
        properties_ = UseMatPsi2(molStr, basisSetAO, basisSetAOandAMBO);
        
    end
    
end
