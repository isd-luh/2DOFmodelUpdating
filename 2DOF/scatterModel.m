classdef scatterModel < handle
    
    properties

        % Lower and upper bounds
        afLB;
        afUB;

        % Output meta model 
        cMetaOutput;

    end
    
    methods(Access = public)
        
        function this = scatterModel(mfInputSamples, mfOutput, afLB, afUB)
            
            this.afLB = afLB;
            this.afUB = afUB;

            % Deduplicate input samples
            [mfuSamples, ind, ~] = unique(mfInputSamples, 'rows');
            
            % Norm input samples
            mfNormeduSamples = (mfuSamples - this.afLB) ./ (this.afUB - this.afLB);
            
            % Interpolation
            this.cMetaOutput = {};
            for iOutput = 1:size(mfOutput, 2)
                this.cMetaOutput{iOutput} = scatteredInterpolant(mfNormeduSamples, mfOutput(ind,iOutput), 'natural');
            end

        end
        
        function mfMetaOutput = predict(this, mfInputSamples)
            
            % Norm input samples
            mfNormedSamples = (mfInputSamples - this.afLB) ./ (this.afUB - this.afLB);
            for iOutput = 1:numel(this.cMetaOutput)
                    mfMetaOutput(:,iOutput) = this.cMetaOutput{iOutput}(mfNormedSamples);
            end

        end
    end
end



