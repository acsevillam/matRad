function p = helper_mvarGauss(model)
%helper_mvarGauss Computes multivariate Gaussian probability for scenario
%model
%   Can be used to test correct scenario probabilities. Also considers used
%   ct phases for probability computation (uncorrelated)
    p = matRad_computeScenarioProbabilities(model.scenarioComponents, ...
        model.scenForProb(:,2:end),model.ctScenProb,model.scenForProb(:,1));
end
