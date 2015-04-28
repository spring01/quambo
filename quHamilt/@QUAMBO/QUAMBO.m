classdef QUAMBO < handle
    
    properties (SetAccess = private)
        
        AOtoQUAMBO;
        
    end
    
    methods
        
        function obj = QUAMBO(input)
            % input section
            func2center = input.funcToCenterAMBO;
            AOtoMO = input.AOtoMO;
            AOandAMBO = input.AOandAMBO;
            numElectrons = input.numElectrons;
            
            % calculate quambo
            MOtoQUAMBO = obj.FindMOtoQUAMBO(AOtoMO, AOandAMBO, numElectrons);
            MOtoQUAMBO = obj.SymmOrthoNShells(MOtoQUAMBO, func2center);
            AOtoQUAMBO_ = AOtoMO * MOtoQUAMBO;                              % <AO\MO><MO\QUAMBO> = <AO\QUAMBO>
            
            % record AOtoQUAMBO
            obj.AOtoQUAMBO = AOtoQUAMBO_;
        end
        
        function newRHF = NewRHF(obj, rhf)
            trans = obj.AOtoQUAMBO;
            prop.overlapMat = trans' * rhf.overlapMat * trans;
            prop.kineticMat = trans' * rhf.kineticMat * trans;
            prop.corePotentialMats = zeros([size(prop.overlapMat), size(rhf.corePotentialMats,3)]);
            for iatom = 1:size(prop.corePotentialMats,3)
                prop.corePotentialMats(:,:,iatom) = trans' * rhf.corePotentialMats(:,:,iatom) * trans;
            end
            prop.twoElecIntegrals = TransformTensor4(rhf.twoElecIntegrals, trans);
            prop.nucRepEnergy = rhf.nucRepEnergy;
            prop.numElectrons = rhf.numElectrons;
            newRHF = RHF(prop);
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
        
        function MOtoQUAMBO = SymmOrthoNShells(~, MOtoQUAMBO, func2center)
            tempOverlap = MOtoQUAMBO' * MOtoQUAMBO;
            for center = unique(func2center)
                funcsOnCenter = find(func2center==center);
                if(length(funcsOnCenter) > 1) % orthogonalize 2s+2p
                    funcs = funcsOnCenter(2:5);
                    MOtoQUAMBO(:,funcs) = MOtoQUAMBO(:,funcs) ...
                        / sqrtm(tempOverlap(funcs,funcs));
                end
                if(length(funcsOnCenter) > 5) % orthogonalize 3s+3p
                    funcs = funcsOnCenter(6:9);
                    MOtoQUAMBO(:,funcs) = MOtoQUAMBO(:,funcs) ...
                        / sqrtm(tempOverlap(funcs,funcs));
                end
                if(length(funcsOnCenter) > 9)
                    throw(MException('QUAMBO:SymmOrthoNShells','Cannot cope with 4th row (and down) for now.'));
                end
            end
        end

    end
    
    methods (Static)
        
        [quambo, matpsi2AO] = MatPsi2Interface(molCart, basisSetInfo);
        
    end
    
end
