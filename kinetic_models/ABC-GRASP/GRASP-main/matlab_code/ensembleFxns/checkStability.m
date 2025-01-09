function isModelValid = checkStability(ensemble,models,strucIdx,eigThreshold)
% Given a model, checks if the largest real part of the jacobian
% eigenvalues is greater than the given threshold (*eigThreshold*), if so 
% the model is considered invalid.
%
%
% USAGE:
%
%    isModelValid = checkStability(ensemble, models, strucIdx, eigThreshold)
%
% INPUT:
%    ensemble (struct):       model ensemble, see buildEnsemble for fields description
%    models (struct):         model, see initialSampler for fields description
%    strucIdx (int):          ID of the model structure considered
%    eigThreshold (double):	  threshold for positive eigenvalues' real part
%
% OUTPUT:
%    isModelValid (logical):	whether or not the model is valid
%
% .. Authors:
%       - Pedro Saa     2018 original code
%       - Marta Matos	2019 refactored code to be used in initialSampler

% Optimization & simulation parameters
kineticFxn   = str2func(ensemble.kineticFxn{strucIdx});
freeVars     = numel(ensemble.freeVars);
Sred         = ensemble.Sred;
numFluxes    = numel(ensemble.fluxRef);
ix_mets      = 1:numel(ensemble.metsBalanced);
ix_enz       = ix_mets(end)+1:freeVars;
xconst       = ones(numel(ensemble.metsFixed), numel(ix_mets));

% fprintf("kineticFxn\n")
% disp(kineticFxn)
% fprintf("freeVars\n")
% disp(freeVars)
% fprintf("Sred\n")
% disp(Sred)
% fprintf("numFluxes\n")
% disp(numFluxes)
% fprintf("ix_mets\n")
% disp(ix_mets)
% fprintf("ix_enz\n")
% disp(ix_enz)
% fprintf("xconst\n")
% disp(xconst)

% Main loop
hstep = 1e-10;              % Step size for control coefficient computations
xopt = ones(freeVars,1);

% fprintf("ones(freeVars,1)\n")
% disp(ones(freeVars,1))

% Define reference state
xref = xopt(ix_mets);
Eref = xopt(ix_enz);

% fprintf("xref\n")
% disp(xref)
% fprintf("Eref\n")
% disp(Eref)

% Define step length to perturb metabolite concentrations
hstep_x = hstep*xref;
xmets   = repmat(xref,1,numel(xref)) + 1i*diag(hstep_x);
xenz    = repmat(Eref,1,numel(xref));
xstep   = [xmets;xenz];

% fprintf("hstep_x\n")
% disp(hstep_x)
% fprintf("xmets\n")
% disp(xmets)
% fprintf("xenz\n")
% disp(xenz)
% fprintf("xstep\n")
% disp(xstep)
% 
% 
% fprintf("kineticFxn")
% disp(kineticFxn)
% save(fullfile('saved', 'kineticFxn.mat'),"kineticFxn")
% 
% fprintf("xstep")
% disp(xstep)
% save(fullfile('saved', 'xstep.mat'),"xstep")
% 
% xstepimag = imag(xstep);
% fprintf("xstepimag")
% disp(xstepimag)
% 
% fprintf("xconst")
% disp(xconst)
% save(fullfile('saved', 'xconst.mat'),"xconst")
% 
% fprintf("models")
% disp(models)
% save(fullfile('saved', 'models.mat'),"models")
% 
% ensemble_fixedExch = ensemble.fixedExch(:,1);
% fprintf("ensemble_fixedExch")
% disp(ensemble_fixedExch)
% save(fullfile('saved', 'ensemble_fixedExch.mat'),"ensemble_fixedExch")
% 
% fprintf("Sred")
% disp(Sred)
% save(fullfile('saved', 'Sred.mat'),"Sred")
% 
% ensemble_kinInactRxns = ensemble.kinInactRxns;
% fprintf("ensemble_kinInactRxns")
% disp(ensemble_kinInactRxns)
% save(fullfile('saved', 'ensemble_kinInactRxns.mat'),"ensemble_kinInactRxns")
% 
% ensemble_subunits = ensemble.subunits{strucIdx};
% fprintf("ensemble_subunits")
% disp(ensemble_subunits)
% save(fullfile('saved', 'ensemble_subunits.mat'),"ensemble_subunits")
% 
% flag = 0;
% fprintf("flag")
% disp(flag)
% save(fullfile('saved', 'flag.mat'),"flag")


% Simulate flux for metabolite perturbation
simFlux = feval(kineticFxn,xstep,xconst,models,ensemble.fixedExch(:,1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);

% fprintf('simFlux:\n')
% disp(simFlux)
% save("simulatedFlux.mat","simFlux")

% Compute elasticiy matrices
E_x_abs  = -(imag(simFlux')./hstep_x(:,ones(1,numFluxes)))';
% fprintf('E_x_abs\n')
% disp(E_x_abs)
% save("E_x_abs.mat","E_x_abs")

% Compute Jacobian eigenvalues
jacobian   = Sred(ensemble.metsLI,:)*E_x_abs(:,ensemble.metsLI);

% fprintf('ensemble.metsLI\n')
% disp(ensemble.metsLI)
% fprintf('Sred\n')
% disp(Sred)
% save("Sred.mat","Sred")
% fprintf('Sred(ensemble.metsLI,:)\n')
% disp(Sred(ensemble.metsLI,:))
% fprintf('jacobian\n')
% disp(jacobian)
% save("jacobian.mat","jacobian")

try
    eigenvalues = eig(jacobian);
catch ME
    if strcmp(ME.identifier, 'MATLAB:eig:inputMustBeSquareStandard')
        error('The jacobian matrix is not square. Hint: check if the correct metabolites are set as constants.');
    else
        rethrow(ME)
    end
end

% fprintf('Eigenvalues:\n')
% disp(eigenvalues)

% Look for positive real eigenvalues
maxEigenvalue = max(real(eigenvalues));

isModelValid = true;

if maxEigenvalue > eigThreshold
    isModelValid = false;
end

end





