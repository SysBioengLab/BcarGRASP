function simulateFluxes(ensembleDir,strucIdx,dataDir,conditionsN,parallelSampling,coresN,solver,backupDir)

%% Load models
ensembleFile = load(ensembleDir);
ensemble = ensembleFile.ensembleFiltered;
addKineticFxnsToPath(ensemble)

%% Load enzyme and metabolites values for simulation

protData = readtable(dataDir,'Sheet','protData');
metsData = readtable(dataDir,'Sheet','metsData');

protOrder = zeros(size(protData,1),1);
metsOrder = zeros(size(metsData,1),1);

for ix = 1:length(protOrder)
    protOrder(ix,1) = find(strcmp(ensemble.rxns,['r_',protData{ix,1}{:}]));
end
for ix = 1:length(metsOrder)
    metsOrder(ix,1) = find(strcmp(ensemble.mets,['m_',metsData{ix,1}{:}]));
end

tempProtData = protData;
tempMetsData = metsData;
for ix = 1:length(protOrder)
    tempProtData(ix,:) = protData(protOrder(ix,1),:);
end
for ix = 1:length(metsOrder)
    tempMetsData(ix,:) = metsData(metsOrder(ix,1),:);
end

protData = tempProtData;
metsData = tempMetsData;

protDataMin = zeros(length(protOrder),conditionsN);
protDataMean = zeros(length(protOrder),conditionsN);
protDataMax = zeros(length(protOrder),conditionsN);
metsDataMin = zeros(length(metsOrder),conditionsN);
metsDataMean = zeros(length(metsOrder),conditionsN);
metsDataMax = zeros(length(metsOrder),conditionsN);

for ix = 1:conditionsN
    protDataMin(:,ix) = protData{:,3*(ix)-1};
    protDataMean(:,ix) = protData{:,3*(ix)};
    protDataMax(:,ix) = protData{:,3*(ix)+1};
    metsDataMin(:,ix) = metsData{:,3*(ix)-1};
    metsDataMean(:,ix) = metsData{:,3*(ix)};
    metsDataMax(:,ix) = metsData{:,3*(ix)+1};
end


%% Initialize parameters

% Solver parameters (FMINCON)
options = optimset('Display','off','Algorithm','sqp','MaxIter',1e4,'TolFun',1e-11,'TolX',1e-10);


%% Perform flux simulation of all the specified conditions

xconst = ones(size(ensemble.metsFixed,1), 1);
kineticFxn = str2func(ensemble.kineticFxn{strucIdx});

simFluxes = cell(conditionsN,1);
xopts = cell(conditionsN,1);
fmins = cell(conditionsN,1);
times = cell(conditionsN,1);
stabilities = cell(conditionsN,1);

% Check if there are pool constraints
if ~isempty(ensemble.poolConst)
    for ix = 1:numel(ensemble.poolConst)
        A{ix} = ensemble.poolConst{ix}(1:numel(ensemble.metsActive));      % extract rhs of from pool constraint matrix
        b{ix} = ensemble.poolConst{ix}(end);
    end
else
    A = [];                                                                % inequality constraints matrix
    b = [];                                                                % rhs ineequality constraints
end
Aeq    = [];                                                               % equality constraints matrix
beq    = [];                                                               % rhs equality constraints
x0     = [metsDataMean;protDataMean; ones(numel(ensemble.kinInactRxns), conditionsN)];                    % initial guess
lb     = [metsDataMin;protDataMin; ones(numel(ensemble.kinInactRxns), conditionsN)];                      % lower bounds on free vars
ub     = [metsDataMax;protDataMax; ones(numel(ensemble.kinInactRxns), conditionsN)];                      % upper bounds on free vars
nlcons = [];
ixTemp = 1;

if isfile(backupDir)
    load(backupDir,"-regexp","^(?!coresN|parallelSampling|solver)...")
end
if ixTemp > conditionsN
    return
end
for ix = ixTemp:conditionsN
    disp(['Simulating Condition: ',num2str(ix)])
    if parallelSampling
        xoptsAux = zeros(size(ensemble.populations.models,2),size(x0,1));
        fminsAux = zeros(size(ensemble.populations.models,2),1);
        timesAux = zeros(size(ensemble.populations.models,2),1);
        simFluxesAux = zeros(size(ensemble.populations.models,2),length(ensemble.rxns));
        stabilitiesAux = zeros(size(ensemble.populations.models,2),1);
        parpool(coresN);

        ensembleFreeVars = ensemble.freeVars;
        ensemblePopulationsModels = ensemble.populations.models;
        ensembleFixedExch= ensemble.fixedExch(:,1);
        ensembleSred = ensemble.Sred;
        ensembleKinInactRxns = ensemble.kinInactRxns;
        ensembleSubunits = ensemble.subunits{strucIdx};
        ensembleDescription = ensemble.description;

        lbTemp = lb(:,ix);
        ubTemp = ub(:,ix);
        x0Temp = x0(:,ix);

        parfor jx = 1:size(ensemble.populations.models,2)
            disp(['Simulating model: ',num2str(jx)])
            % NLOPT call
            if strcmp(solver,'NLOPT')
                opt = struct;
                opt.algorithm = 40; 									   			   % (NLOPT_LD_MMA, descartado) 11(NLOPT_LD_LBFGS), options: 40(NLOPT_LD_SLSQP), 13(NLOPT_LD_VAR1), 14(NLOPT_LD_VAR2)
                opt.ftol_abs  = 1e-11;
                opt.xtol_abs  = 1e-10*ones(1,numel(ensembleFreeVars));
                opt.maxeval   = 1e5;
                opt.maxtime   = 7200;
                opt.min_objective = @(x) feval(kineticFxn,x,xconst,ensemblePopulationsModels(jx),ensembleFixedExch,ensembleSred,ensembleKinInactRxns,ensembleSubunits,1); %does not support fixed exchanges yet
                opt.lower_bounds  = lbTemp;
                opt.upper_bounds  = ubTemp;
                if contains(ensembleDescription,'b_car_rejection_model_detailed') %% ONLY FOR DETAILED MODELS
                    opt.fc{1,1} = (@(x) detailedEnzymeConstrained(x,38,39,1));
                    opt.fc{1,2} = (@(x) detailedEnzymeConstrained(x,38,39,-1));
                    opt.fc{1,3} = (@(x) detailedEnzymeConstrained(x,42,43,1));
                    opt.fc{1,4} = (@(x) detailedEnzymeConstrained(x,42,43,-1));
                    opt.fc{1,5} = (@(x) detailedEnzymeConstrained(x,44,45,1));
                    opt.fc{1,6} = (@(x) detailedEnzymeConstrained(x,44,45,-1));
                    opt.fc{1,7} = (@(x) detailedEnzymeConstrained(x,46,47,1));
                    opt.fc{1,8} = (@(x) detailedEnzymeConstrained(x,46,47,-1));
                    opt.fc{1,9} = (@(x) detailedEnzymeConstrained(x,46,48,1));
                    opt.fc{1,10} = (@(x) detailedEnzymeConstrained(x,46,48,-1));
                    opt.fc{1,11} = (@(x) detailedEnzymeConstrained(x,46,49,1));
                    opt.fc{1,12} = (@(x) detailedEnzymeConstrained(x,46,49,-1));
                    opt.fc_tol = 1e-6*ones(1,12);
                end
                tic
                [xopt,fmin, ~] = nlopt_optimize(opt,x0Temp);
                time = toc;
                xoptsAux(jx,:) = xopt';
                fminsAux(jx,:) = fmin;
                timesAux(jx,:) = time;
            % FMINCON call
            else
                Aeq_fmin = [];
                beq_fmin = [];
                if contains(ensembleDescription,'b_car_rejection_model_detailed') %% ONLY FOR DETAILED MODELS
                    Aeq_fmin = zeros(6,51);
                    beq_fmin = zeros(6,1);
                    Aeq_fmin(1,38) = 1;
                    Aeq_fmin(1,39) = -1;
                    Aeq_fmin(2,42) = 1;
                    Aeq_fmin(2,43) = -1;
                    Aeq_fmin(3,44) = 1;
                    Aeq_fmin(3,45) = -1;
                    Aeq_fmin(4,46) = 1;
                    Aeq_fmin(4,47) = -1;
                    Aeq_fmin(5,46) = 1;
                    Aeq_fmin(5,48) = -1;
                    Aeq_fmin(6,46) = 1;
                    Aeq_fmin(6,49) = -1;
                end
                tic;
                [xopt,fmin,~] = fmincon(kineticFxn,x0Temp,[],[],Aeq_fmin,beq_fmin,lbTemp,ubTemp,[],options,xconst,ensemblePopulationsModels(jx),ensembleFixedExch,ensembleSred,ensembleKinInactRxns,ensembleSubunits,1);
                time = toc;
                xoptsAux(jx,:) = xopt';
                fminsAux(jx,:) = fmin;
                timesAux(jx,:) = time;
            end
            simulatedFlux = feval(kineticFxn,xopt,xconst,ensemblePopulationsModels(jx),ensembleFixedExch,ensembleSred,ensembleKinInactRxns,ensembleSubunits,0);
            simFluxesAux(jx,:) = simulatedFlux';
            stabilitiesAux(jx,:) = checkStabilityRejection(ensemble,ensemble.populations.models(jx),strucIdx, ensemble.eigThreshold,xopt);
        end
        xopts{ix,1} = xoptsAux;
        fmins{ix,1} = fminsAux;
        times{ix,1} = timesAux;
        simFluxes{ix,1} = simFluxesAux;
        stabilities{ix,1} = stabilitiesAux;
        delete(gcp)
        
    else
        for jx = 1:size(ensemble.populations.models,2)
            disp(['Simulating model: ',num2str(jx)])
            % NLOPT call
            if strcmp(solver,'NLOPT')
                opt = struct;
                opt.algorithm = 40; 									   			   % (NLOPT_LD_MMA, descartado) 11(NLOPT_LD_LBFGS), options: 40(NLOPT_LD_SLSQP), 13(NLOPT_LD_VAR1), 14(NLOPT_LD_VAR2)
                opt.ftol_abs  = 1e-11;
                opt.xtol_abs  = 1e-10*ones(1,numel(ensemble.freeVars));
                opt.maxeval   = 1e5;
                opt.maxtime   = 7200;
                opt.min_objective = @(x) kineticFxn(x,xconst,ensemble.populations.models(jx),ensemble.fixedExch(:,1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},1); %does not support fixed exchanges yet
                opt.lower_bounds  = lb(:,ix);
                opt.upper_bounds  = ub(:,ix);
                if contains(ensemble.description,'b_car_rejection_model_detailed') %% ONLY FOR DETAILED MODELS
                    opt.fc{1,1} = (@(x) detailedEnzymeConstrained(x,38,39,1));
                    opt.fc{1,2} = (@(x) detailedEnzymeConstrained(x,38,39,-1));
                    opt.fc{1,3} = (@(x) detailedEnzymeConstrained(x,42,43,1));
                    opt.fc{1,4} = (@(x) detailedEnzymeConstrained(x,42,43,-1));
                    opt.fc{1,5} = (@(x) detailedEnzymeConstrained(x,44,45,1));
                    opt.fc{1,6} = (@(x) detailedEnzymeConstrained(x,44,45,-1));
                    opt.fc{1,7} = (@(x) detailedEnzymeConstrained(x,46,47,1));
                    opt.fc{1,8} = (@(x) detailedEnzymeConstrained(x,46,47,-1));
                    opt.fc{1,9} = (@(x) detailedEnzymeConstrained(x,46,48,1));
                    opt.fc{1,10} = (@(x) detailedEnzymeConstrained(x,46,48,-1));
                    opt.fc{1,11} = (@(x) detailedEnzymeConstrained(x,46,49,1));
                    opt.fc{1,12} = (@(x) detailedEnzymeConstrained(x,46,49,-1));
                    opt.fc_tol = 1e-6*ones(1,12);
                end
                tic
                [xopt,fmin, retcode] = nlopt_optimize(opt,x0(:,ix));
                time = toc;
                xopts{ix,1}(jx,:) = xopt';
                fmins{ix,1}(jx,:) = fmin;
                times{ix,1}(jx,:) = time;
            % FMINCON call
            else
                Aeq_fmin = [];
                beq_fmin = [];
                if contains(ensemble.description,'b_car_rejection_model_detailed') %% ONLY FOR DETAILED MODELS
                    Aeq_fmin = zeros(6,51);
                    beq_fmin = zeros(6,1);
                    Aeq_fmin(1,38) = 1;
                    Aeq_fmin(1,39) = -1;
                    Aeq_fmin(2,42) = 1;
                    Aeq_fmin(2,43) = -1;
                    Aeq_fmin(3,44) = 1;
                    Aeq_fmin(3,45) = -1;
                    Aeq_fmin(4,46) = 1;
                    Aeq_fmin(4,47) = -1;
                    Aeq_fmin(5,46) = 1;
                    Aeq_fmin(5,48) = -1;
                    Aeq_fmin(6,46) = 1;
                    Aeq_fmin(6,49) = -1;
                end
                tic;
                [xopt,fmin,retcode] = fmincon(kineticFxn,x0(:,ix),[],[],Aeq_fmin,beq_fmin,lb(:,ix),ub(:,ix),[],options,xconst,ensemble.populations.models(jx),ensemble.fixedExch(:,1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},1);
                time = toc;
                xopts{ix,1}(jx,:) = xopt';
                fmins{ix,1}(jx,:) = fmin;
                times{ix,1}(jx,:) = time;
            end
            simulatedFlux = feval(kineticFxn,xopt,xconst,ensemble.populations.models(jx),ensemble.fixedExch(:,1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);
            simFluxes{ix,1}(jx,:) = simulatedFlux';
            stabilities{ix,1}(jx,:) = checkStabilityRejection(ensemble,ensemble.populations.models(jx),strucIdx, ensemble.eigThreshold,xopt);
        end
    end
    ixTemp = ix+1;
    save(backupDir) %
end
