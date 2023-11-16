classdef class2DOF < handle

    properties

        % Queue for saving design variables & eigenvalues
        mfQueueResults;

        % Correct eigenvalues
        afEigenvalues

        % Standard distribution eigenvalues (assumed noise)
        afSigmaNoise

        % Method (choices are 'TMCMC', 'MCGO', 'meta-MCGO'
        strMethod

        % Output meta model (eigenvalues)
        metaOutput


    end

    methods

        %% Constructor of class
        function this = class2DOF(afEigenvalues, afSigmaNoise, mfQueueResults, strMethod, metaOutput)

            this.afEigenvalues = afEigenvalues;
            this.afSigmaNoise = afSigmaNoise;
            this.mfQueueResults  = mfQueueResults;
            this.strMethod = strMethod;
            this.metaOutput = metaOutput;

        end

        %% Calculate objective function value
        function cResult = calcObjValue(this, x, ~)

            if strcmp(this.strMethod, 'MCGO') || strcmp(this.strMethod, 'meta-MCGO')
                if isempty(x)
                    cResult = 1;
                    return
                end
            end

            % Calculate current eigenvalues
            if strcmp(this.strMethod, 'MCGO') || strcmp(this.strMethod, 'TMCMC')
                afCurrentEigenvalues = calcEigenvalues(x);
            elseif strcmp(this.strMethod, 'meta-MCGO')
                afCurrentEigenvalues = this.metaOutput(x);
            end
            nDV = numel(this.afEigenvalues);

            % Send design variables and currrent eigenvalues to queue
            send(this.mfQueueResults, [x, afCurrentEigenvalues]);

            % Calculate objective function
            if strcmp(this.strMethod, 'TMCMC')
                for iDV = 1:nDV
                    logl_i(iDV) = (1/this.afSigmaNoise(iDV))^2 * ...
                        (this.afEigenvalues(iDV) - afCurrentEigenvalues(iDV))^2;
                end
                res =  - 0.5 * sum(logl_i);
                cResult = {res, [], []};

            elseif strcmp(this.strMethod, 'MCGO') || strcmp(this.strMethod, 'meta-MCGO')
                for iDV = 1:nDV
                    diffEigenvalues(iDV) = (afCurrentEigenvalues(iDV) - ...
                        this.afEigenvalues(:,iDV))^2;
                end
                cResult = sum(diffEigenvalues);

            end

        end

    end
end