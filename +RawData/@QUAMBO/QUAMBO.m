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
            MOtoQUAMBO = obj.CalculateMOtoQUAMBO(properties.AOtoMO, ...
                properties.AOandAMBO, properties.numElectrons);
            obj.AOtoQUAMBO = properties.AOtoMO * MOtoQUAMBO;                % <AO\MO><MO\QUAMBO> = <AO\QUAMBO>
            
            % integrals
            obj.overlapQUAMBO = MOtoQUAMBO' * MOtoQUAMBO;
            obj.kineticQUAMBO = ...
                obj.AOtoQUAMBO' * properties.kineticAO * obj.AOtoQUAMBO;
            obj.potentialEachCoreQUAMBO = zeros( ...
                [size(obj.kineticQUAMBO), ...
                size(properties.potentialEachCoreAO, 3)]);
            for iAtom = 1:size(properties.potentialEachCoreAO, 3)
                obj.potentialEachCoreQUAMBO(:,:,iAtom) ...
                    = obj.AOtoQUAMBO' ...
                    * properties.potentialEachCoreAO(:,:,iAtom) ...
                    * obj.AOtoQUAMBO;
            end
            obj.twoElecIntegralsQUAMBO = obj.TransformTensor4( ...
                properties.twoElecIntegralsAO, obj.AOtoQUAMBO);
        end
        
    end
    
    methods (Access = private)
        
        function MOtoQUAMBO = CalculateMOtoQUAMBO(~, AOtoMO, AOandAMBO, numElectrons)
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
        
        properties_ = UseMatPsi2(molStr, basisSetAO, basisSetAOandAMBO, packageLocation);
        
    end
    
end
