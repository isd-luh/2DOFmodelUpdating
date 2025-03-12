classdef class2DOF < handle

    properties

        % Queue for saving design variables & eigenvalues
        queueResults;

        % Mean eigenvalues
        meanEigenvalues

        % Standard distribution eigenvalues 
        sigmaEigenvalues

        % Method (choices are 'BMU', 'RDMU', 'MetaRDMU'
        strMethod

        % Output meta model (eigenvalues)
        outputMetaModel


    end

    methods

        %% Constructor of class
        function this = class2DOF(meanEigenvalues, sigmaEigenvalues, queueResults, strMethod, outputMetaModel)

            this.meanEigenvalues = meanEigenvalues;
            this.sigmaEigenvalues = sigmaEigenvalues;
            this.queueResults  = queueResults;
            this.strMethod = strMethod;
            this.outputMetaModel = outputMetaModel;

        end

        %% Calculate objective function value
        function result = calcObjValue(this, x, ~)

            if strcmp(this.strMethod, 'RDMU') || strcmp(this.strMethod, 'MetaRDMU')
                if isempty(x)
                    result = 1;
                    return
                end
            end

            % Calculate current eigenvalues
            if strcmp(this.strMethod, 'RDMU') || strcmp(this.strMethod, 'BMU')
                currentEigenvalues = calcEigenvalues(x);
            elseif strcmp(this.strMethod, 'MetaRDMU')
                currentEigenvalues = this.outputMetaModel(x);
            end
            nDV = numel(this.meanEigenvalues);

            % Send design variables and currrent eigenvalues to queue (only for RDMU approach)
            if ~isempty(this.queueResults)
                send(this.queueResults, [x, currentEigenvalues]);
            end

            % Calculate objective/likelihood function
            if strcmp(this.strMethod, 'BMU')
                for iDV = 1:nDV
                    logl_i(iDV) = (1/this.sigmaEigenvalues(iDV))^2 * ...
                        (this.meanEigenvalues(iDV) - currentEigenvalues(iDV))^2;
                end
                res =  - 0.5 * sum(logl_i);
                result = {res, [], []};

            elseif strcmp(this.strMethod, 'RDMU') || strcmp(this.strMethod, 'MetaRDMU')
                for iDV = 1:nDV
                    diffEigenvalues(iDV) = (currentEigenvalues(iDV) - ...
                        this.meanEigenvalues(:,iDV))^2;
                end
                result = sum(diffEigenvalues);

            end

        end

    end
end