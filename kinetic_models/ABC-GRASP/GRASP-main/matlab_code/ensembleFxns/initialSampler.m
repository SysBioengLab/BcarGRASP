function [isModelValid,models,strucIdx,xopt,tolScore,simulatedFlux,fmins,times] = initialSampler(ensemble, modelI)
% Samples initial ensemble of kinetic models.
%
%
% Checks if sampled models are valid and only returns valid ones.
% A model is considered valid if
%
%  - for all reactions the fluxes and respective Gibbs energies are 
%    compatible;
%  - the real part of the jacobian eigenvalues is lower than the defined 
%    threshold;
%  - the difference between the predicted flux and the reference flux
%    is negligible.
%  - if the rejection sampler is used, the difference between the model 
%    and the data is smaller than the defined tolerance.
%
%
% USAGE:
%
%    [isModelValid, models, strucIdx, xopt, tolScore, simulatedFlux] = initialSampler(ensemble)
%
% INPUT:
%    ensemble (struct):  initialized model ensemble, see initializeEnsemble for a list of all fields in the ensemble struct.
%
% OUTPUT:
%    isModelValid (logical):   whether or not model is valid.
%    models (struct):          sampled model
%
%               * poolFactor (*double vector*) : relative pool proportions of specific metabolites at reference state (e.g., atp/adp)
%               * refFlux (*double vector*)    : reference reaction fluxes (mean)
%               * metConcRef (*double vector*) : reference metabolite concentrations (sampled within thermodynamically feasible ranges)
%               * gibbsTemp (*double vector*)  : Gibbs reactions energies
%               * rxnParams (*struct*)         : reactions parameters
%
%                       * reversibilities (*double vector*)  : sampled elementary reaction reversibilities
%                       * enzymeAbundances (*double vector*) : sampled enzyme intermediates abundances
%                       * modifierElemFlux (*double vector*) : [TODO Pedro]
%                       * branchFactor (*double vector*)     : [TODO Pedro]
%                       * kineticParams (*double vector*)    : reaction kinetic parameters
%    strucIdx (int):                   model structure ID
%    xopt (*double vector*):           predicted metabolite and enzyme concentrations of an experimentally consistent model
%    tolScore (*double vector*):       discrepancy tolerance between predicted and experimental fluxes
%    simulatedFlux (*double vector*):  predicted fluxes
%
% .. Authors:
%       - Pedro Saa         2016 original code
%       - Marta Matos       2018, 2019 generalized it for promiscuous  
%                           reactions and random mechanisms, added model 
%                           validity checks
%       - Nicholas Cowie	2019 added extreme pathways and random flux 
%                           distribution for isoenzymes

%% Initialze parameters
RT       = 8.314*298.15/1e3;                                               % gas constant times the absolute temperature (298.15 K)
massTol  = size(ensemble.Sred,1)*1e-10;								       % #balances*tol^2

% Just so the tests don't crash because these variables were not assigned
xopt = 0;
tolScore = 0;
simulatedFlux = 0;
fmins = 0;
times = 0;

% Figure out NLP solver
if strcmpi(ensemble.solver,'NLOPT')											   % Solver parameters for NLOPT
    opt.algorithm = 40; 									   			   % (NLOPT_LD_MMA, descartado) 11(NLOPT_LD_LBFGS), options: 40(NLOPT_LD_SLSQP), 13(NLOPT_LD_VAR1), 14(NLOPT_LD_VAR2)
    opt.ftol_abs  = 1e-11;
    opt.xtol_abs  = 1e-10*ones(1,numel(ensemble.freeVars));
    opt.maxeval   = 1e5;
    opt.maxtime   = 7200;
    % opt.verbose   = 1;
elseif strcmpi(ensemble.solver,'OPTI_IPOPT')
    opts = optiset('solver','ipopt','display','off','maxiter',1e4,'maxfeval',10000,'maxtime',1000,'maxnodes',10000,'tolrfun',1e-11,'tolafun',1e-11,'tolint',1e-10);
elseif strcmpi(ensemble.solver,'OPTI_NLOPT')
    if strcmp(ensemble.description,'b_car_rejection_model_detailed') %
        opt_time = 9800; %
    else %
        opt_time = 360; %
    end %
    nopts = nloptset('algorithm','LD_SLSQP');
    opts = optiset('solver','nlopt','display','off','maxiter',1e4,'maxfeval',10000,'maxtime',opt_time,'maxnodes',10000,'tolrfun',1e-11,'tolafun',1e-11,'tolint',1e-10,'solverOpts',nopts);
elseif strcmpi(ensemble.solver,'FMINCON') || isempty(ensemble.solver)							   % Solver parameters (FMINCON)
    options = optimset('Display','off','Algorithm','sqp','MaxIter',1e4,'TolFun',1e-11,'TolX',1e-10);
end

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
x0     = [ensemble.metsDataMean;ensemble.protDataMean; ones(numel(ensemble.kinInactRxns), ensemble.numConditions)];                    % initial guess
lb     = [ensemble.metsDataMin;ensemble.protDataMin; ones(numel(ensemble.kinInactRxns), ensemble.numConditions)];                      % lower bounds on free vars
ub     = [ensemble.metsDataMax;ensemble.protDataMax; ones(numel(ensemble.kinInactRxns), ensemble.numConditions)];                      % upper bounds on free vars
nlcons = [];                                                               % nonlinear constraints (not used)

%% Execute Rejection-ABC
acceptanceRate = 1;
counter        = 0;

% Loop until the number of valid particles is reached
while true
    isModelValid = true;
    
    % Update attempt counter
    counter = counter+1;

    % Sample model structure
    strucIdx = randi(ensemble.numStruct);

    % Sample pool parameters (if any)
    if ~isempty(ensemble.poolConst)
        poolFactor = [];
        poolFactor{size(ensemble.poolConst{1},1)} = 0;
        for ix = 1:numel(poolFactor)

            % Generate pool factor ~ Dir(alpha) using independent gamma distributions
            alphaPoolFactor = ensemble.populations(1).probParams(strucIdx).alphaPoolFactor{ix};
            poolFactorTemp  = randg(alphaPoolFactor);
            poolFactor      = poolFactorTemp/sum(poolFactorTemp);

            % Update pool constraint matrix accordingly
            for jx = 1:ensemble.numConditions                        
                A_opt{jx} = A{jx};
                A_opt{jx}(A{jx}~=0) = poolFactor;
            end            
            if (ix==1)
                models(1).poolFactor = poolFactor;
            else
                models(1).poolFactor = [models(1).poolFactor;poolFactor];
            end
        end
    else
        models(1).poolFactor = [];
    end

    % Randomly distribute flux between isoenzymes

    if ~isempty(ensemble.uniqueIso)
        for xi = 1:size(ensemble.uniqueIso,1)
            group = find(strcmp(ensemble.isoenzymes,ensemble.uniqueIso{xi}));
            splitFactor = zeros(size(group,1),1);
            totalFlux = sum(ensemble.fluxPoints(group, modelI));
            for yi = 1:size(splitFactor,1)
                splitFactor(yi) = randg();
            end
            splitFactor = splitFactor./sum(splitFactor);
            ensemble.fluxPoints(group, modelI) = splitFactor.*totalFlux;
        end
    end
  
    models(1).refFlux =  ensemble.fluxPoints(:, modelI); 
    models(1).fixedExch = models(1).refFlux(ensemble.kinInactRxns,:);
    assert(all(abs(ensemble.Sred * models.refFlux) <10^-8), "Your model doesn\'t seem to be at steady-state. Sred * fluxRef != 0");

    % Determine gibbs free energy of reaction
    models(1).gibbsTemp = ensemble.gibbsEnergies(:, modelI);
    models(1).metConcRef = ensemble.metConcRef(:, modelI); 

    % Sample Reversibilities
    [ensemble, models, isModelValid] = sampleGeneralReversibilities(ensemble, models, RT, strucIdx);
    if ~isModelValid
        break;
    end

    % Sample enzyme abundances
    models = sampleEnzymeAbundances(ensemble,models,strucIdx);

    % Sample modifier elementary fluxes (positions are given where exp(R)=1)
    models = sampleModifierElemFluxes(ensemble, models, strucIdx);

    % Calculate rate parameters for allosteric reaction part;
    [ensemble, models] = sampleAllostery(ensemble, models, strucIdx);

    %thermoCounter   = 1;
    for activRxnIdx = 1:numel(ensemble.kinActRxns)
        %disp(ensemble.rxns(ensemble.kinActRxns(activRxnIdx)));

        % Case 1: Diffusion and Exchanges
        if strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')||...
                strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange')
            models(1).rxnParams(ensemble.kinActRxns(activRxnIdx)).kineticParams =  models(1).refFlux(ensemble.kinActRxns(activRxnIdx));

        % Case 2: Enzymatic reactions
        else

            % Check whether the reaction is mass action
            if strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')
                reactionFlux = models(1).refFlux(ensemble.kinActRxns(activRxnIdx));
                gibbsTemp = models(1).gibbsTemp(ensemble.kinActRxns(activRxnIdx));
                models(1).rxnParams(activRxnIdx).kineticParams = [1,exp(gibbsTemp/RT)]*reactionFlux/(1-exp(gibbsTemp/RT));
                continue;
            end

            promiscRxnsList = ensemble.promiscuity{strucIdx}{ensemble.kinActRxns(activRxnIdx)};
            reverTemp = ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)};
            reactionFlux = ensemble.reactionFluxAllosteric(ensemble.kinActRxns(activRxnIdx));
            randomEnzymesR = models(1).rxnParams(activRxnIdx).enzymeAbundances';
            extremePathways = ensemble.extremePathways{strucIdx}{activRxnIdx};


            % Sample branching factor (if necessary)
            branchFactor = 1;
            Nelem        = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};

            % If the reaction is promiscuous
            if size(promiscRxnsList,1) > 0
                rxnIsPromiscuous = true;
                branchFactor = ensemble.reactionFluxAllosteric(promiscRxnsList)';
                branchFactor = branchFactor/max(branchFactor);

                % If the promiscuous reactions share common steps
                if sum(sum(Nelem)) > size(Nelem,1)
                    reactionFlux = sum(ensemble.reactionFluxAllosteric(promiscRxnsList));

                % If the promiscuous reactions do not share any common steps
                else
                    reactionFlux = max(ensemble.reactionFluxAllosteric(promiscRxnsList));
                end

                % For non promiscuous reactions
            else
                rxnIsPromiscuous = false;
                if (size(extremePathways,2)>1)
                    branchFactor = zeros(1,size(extremePathways,2));
                    for ix = 1:size(extremePathways,2)
                        aBranch            = randg(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaBranchFactor(ix,:));
                        branchFactor(1, ix) = aBranch;
                    end
                    branchFactor = branchFactor/sum(branchFactor);
                end
            end
            models(1).rxnParams(activRxnIdx).branchFactor = branchFactor';

            % Get modifier elementary fluxes (positions are given were exp(R)=1)
            modifierElemFlux = models(1).rxnParams(activRxnIdx).modiferElemFlux';
            % VI. Calculate rate parameters
            %disp(ensemble.rxns{ensemble.kinActRxns(activRxnIdx),strucIdx});
            forwardFlux    = ensemble.forwardFlux{ensemble.kinActRxns(activRxnIdx),strucIdx};
            models(1).rxnParams(activRxnIdx).kineticParams = ...
                calculateKineticParams(reverTemp,forwardFlux,reactionFlux,randomEnzymesR,extremePathways,branchFactor,modifierElemFlux,rxnIsPromiscuous);
%             fprintf("models(1).rxnParams(activRxnIdx).kineticParams")
%             disp(models(1).rxnParams(activRxnIdx).kineticParams)
        end
    end
    
    % Test model consistency
    xconst = ones(size(ensemble.metsFixed,1), 1);
    kineticFxn = str2func(ensemble.kineticFxn{strucIdx});
    testFlux   = feval(kineticFxn,ones(size(ensemble.freeVars,1),1),xconst,models,models(1).fixedExch,ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);
    
%     fprintf("testFlux\n")
%     disp(testFlux)
%     fprintf("models\n")
%     disp(models)
%     save("models.mat","models")
%     fprintf("models(1).refFlux\n")
%     disp(models(1).refFlux)
%     save("test_ensemble.mat","ensemble")

    % If the model is consistent continue
    if any(abs(testFlux-models(1).refFlux)>1e-6) || any(isnan(testFlux))
        isModelValid = false;
        disp(['There are consistency problems during the reaction sampling. Model ID: ',num2str(modelI)]);
        return
    end
    
    % Test if the real part of the jacobian's eigenvalue is greater than
    % threshold
    isModelValid = checkStability(ensemble,models,strucIdx, ensemble.eigThreshold);
    if ~isModelValid
        disp(['There are eigenvalues larger than ', num2str(ensemble.eigThreshold), '. Model ID: ',num2str(modelI)]);
        return
    end
    

    % Check sampling mode. For the GRASP mode, no need to simulate
    if strcmpi(ensemble.sampler,'GRASP'); break;

        % For the remaining modes, we need to simulate the model in the
        % experimental conditions
    elseif ~strcmpi(ensemble.sampler,'GRASP') && isModelValid

        % Simulate fluxes
        tolScore      = 10000*ones(1, ensemble.numConditions);
        simulatedFlux = zeros(numel(ensemble.activeRxns),ensemble.numConditions);
        xopt          = zeros(size(x0,1),ensemble.numConditions);
        fmins = 1000*ones(1,ensemble.numConditions);
        times = zeros(1, ensemble.numConditions);

        % Solver call (OPTI Toolbox not implemented yet)
        for ix = 1:ensemble.numConditions

            % NLOPT call
            if strcmpi(ensemble.solver,'NLOPT')

                % Define anonymous function with objective and bound
                % constraints
                opt.min_objective = @(x) kineticFxn(x,xconst,models,ensemble.fixedExch(:,ix+1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},1);
                opt.lower_bounds  = lb(:,ix);
                opt.upper_bounds  = ub(:,ix);

                % Solves S*v(k,X) = 0; s.t. A*X <= beq, lb < X <ub, with extra constraints (e.g., pool or ratio constraints). Otherwise solve solve S*v(k,X) = 0; s.t. lb < X <ub, with no extra constraints
                if ~isempty(ensemble.poolConst)
                    for jx = 1:numel(ensemble.poolConst)
                        opt.fc{1,2*jx-1} = (@(x) poolConstraintFxn(x,[A_opt{jx},zeros(1,numel(x0(:,ix))-numel(ensemble.metsActive))],b{jx}));
                        opt.fc{1,2*jx}   = (@(x) poolConstraintFxn(x,[-A_opt{jx},zeros(1,numel(x0(:,ix))-numel(ensemble.metsActive))],-b{jx}));
                    end
                    opt.fc_tol = 1e-6*ones(1,2*numel(ensemble.poolConst));
                end
                % Modification to include enzymes with same proportions
                % ONLY for b_car model
                if contains(ensemble.description,'b_car_rejection_model_detailed')
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
                % End of modifications only for b_car models
                tic
                [xopt(:,ix),fmin, retcode] = nlopt_optimize(opt,x0(:,ix));
                time = toc;
                fmins(:,ix) = fmin;
                times(:,ix) = time;
                if retcode < 0
                    disp(['Model number ',mat2str(modelI),' failed ABC optimization. Error code: ', retcode])
                    isModelValid = false;
                    return
                    % error(['The ABC optimization with NLOPT was not successful. Error code: ', retcode, '. For more information, see "Return values" in https://nlopt.readthedocs.io/en/latest/NLopt_Reference/']);
                end

                % OPTI
            elseif strcmpi(ensemble.solver,'OPTI_IPOPT') || strcmpi(ensemble.solver,'OPTI_NLOPT')
                fun_opti = @(x) kineticFxn(x,xconst,models,ensemble.fixedExch(:,ix+1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},1);
                Aeq_opti = [];
                beq_opti = [];
                x0_opti = x0(:,ix);
                lb_opti = lb(:,ix);
                ub_opti = ub(:,ix);
                if ~isempty(ensemble.poolConst)
                    Aeq_opti = [A_opt{ix},zeros(1,numel(x0(:,ix))-numel(ensemble.metsActive))];
                    beq_opti = b{ix};
                end
                Opt = opti('fun',fun_opti,'eq',Aeq_opti,beq_opti,'bounds',lb_opti,ub_opti,'options',opts);
                tic;
                [xopt(:,ix),fmin,exitflag,info] = solve(Opt,x0_opti);
                time = toc;
                fmins(:,ix) = fmin;
                times(:,ix) = time;
                if exitflag < 0
                    display(info)
                    disp(['Model number ',mat2str(modelI),' failed ABC optimization. Error code: ', exitflag])
                    isModelValid = false;
                    return
                    % error(['The ABC optimization with OPTI was not succesful. Error code: ',exitflag,'. For more information, check the OPTI documentation on exitflags.'])
                end
            else
                % Solves S*v(k,X) = 0; s.t. Aeq*X = beq, lb < X <ub, with extra constraints (e.g., pool or ratio constraints). Otherwise solve solve S*v(k,X) = 0; s.t. lb < X <ub, with no extra constraints
                Aeq_fmin = [];
                beq_fmin = [];
                if ~isempty(ensemble.poolConst)                    
                    Aeq_fmin = [A_opt{ix},zeros(1,numel(x0(:,ix))-numel(ensemble.metsActive))];
                    beq_fmin = b{ix};
                end
                tic;
                [xopt(:,ix),fmin,retcode] = fmincon(kineticFxn,x0(:,ix),[],[],Aeq_fmin,beq_fmin,lb(:,ix),ub(:,ix),[],options,xconst,models,ensemble.fixedExch,ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},1);
                time = toc;
                fmins(:,ix) = fmin;
                times(:,ix) = time;
                if retcode < 0
                    disp(['Model number ',mat2str(modelI),' failed ABC optimization. Error code: ', retcode])
                    isModelValid = false;
                    return
                    % error(['The ABC optimization with fmincon was not successful. Error code: ', retcode, '. For more information, check the Matlab documentation on fmincon.']);
                end
            
            end

            % Check mass balance consistency
%             if (fmin<massTol)

                % Simulate fluxes if the system is mass-balanced
                simulatedFlux(:,ix) = feval(kineticFxn,xopt(:,ix),xconst,models,ensemble.fixedExch(:,ix+1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);
                % Calculate discrepancy score
                tolScore(ix) = max(sqrt(mean(((simulatedFlux(ensemble.freeFluxes,ix)-ensemble.simWeights(:,ix))./ensemble.simWeights(:,ix)).^2)));
%                 disp(['Condition ',num2str(ix),' score: ',num2str(tolScore(ix)), ', fmin: ',num2str(fmin),', time: ',num2str(time)])
                % Check tolerance inmediately for this condition
%                 if (tolScore(ix)>ensemble.tolerance)
%                     disp(['The model score is higher than the defined tolerance: ', num2str(tolScore(ix)), ', tolerance: ', num2str(ensemble.tolerance), '.',' Failed at idx: ',num2str(ix)]);
%                     isModelValid = false;
%                     return;
%                end
%             else
%                 isModelValid = false;
%                 disp(['Condition ',num2str(ix),' fmin: ',num2str(fmin),'> massTol: ',num2str(massTol)])
%                 return;
%             end
        end
        disp(['Model number ',mat2str(modelI),' performance:'])
        disp(['fmin: ',mat2str(fmins),', massTol: ',num2str(massTol)])
        disp(['model scores: ', mat2str(round(tolScore,3)), ', tolerance: ', num2str(ensemble.tolerance)]);
        disp(['times: ',mat2str(times)])
        if any(fmins > massTol)
            %disp(['Conditions fmin: ',mat2str(fmins),'> massTol: ',num2str(massTol)])
            isModelValid = false;
        end

        if any(tolScore > ensemble.tolerance)
            %disp(['The model score is higher than the defined tolerance: ', mat2str(round(tolScore,3)), ', tolerance: ', num2str(ensemble.tolerance)]);
            isModelValid = false;
        end

        if isModelValid
            ensemble_test = ensemble;
            ensemble_test.populations(1).strucIdx(1)  = strucIdx;                                                                           % model structures
            ensemble_test.populations(1).tolScore(1,:)  = tolScore;                                                                           % tolerance score
            ensemble_test.populations(1).xopt{1}      = xopt;                                                                               % optimal value found
            ensemble_test.populations(1).simFluxes{1} = simulatedFlux;                                                                          % simulated fluxes
            ensemble_test.populations(1).models    = models; 
            saveMCAMatrices = 1;
            mcaResults = controlAndResponseAnalysis(ensemble_test,saveMCAMatrices);
            for jx = 1:(ensemble.numConditions)
                isModelValid = checkStabilityRejection(ensemble,models,strucIdx, ensemble.eigThreshold,xopt(:,jx));
                if ~isModelValid
                    disp(['Condition ',num2str(jx),' has failed Stability Check of the defined threshold'])
                end
            end
            for jx = 1:(ensemble.numConditions+1)
                if isempty(mcaResults.xControl{jx}) || isempty(mcaResults.vControl{jx})
                    disp(['Condition ',num2str(jx-1),' failed MCA Analysis. Model is rejected'])
                    isModelValid = false;
                    return
                end
            end
        end

        if ~isModelValid
            return
        end


        % Compute tolerance, acceptance rate and break
        if isModelValid
            acceptanceRate = 1/counter;
            return;

            % Delete model if not accurate or mass-balanced
        else
            models = [];
        end
    end
end

% Save results and write progress to a temp file (except for the GRASP mode)
% if ~strcmpi(ensemble.sampler,'GRASP')
%     try
%         load progress.txt
%         progress = [progress(1)+1;progress(4)/(progress(1)+1);progress(5)/(progress(1)+1);progress(4)+acceptanceRate;progress(5)+tolScore];
%         save progress.txt -ascii progress
%         save(['temp/particle_',num2str(progress(1)),'.mat'],'models','strucIdx','xopt','tolScore','simulatedFlux');
% 
%         % If another worker is writing on the file, wait a brief random time
%     catch
%         pause(randi(2)*rand(1));
%         load progress.txt
%         progress = [progress(1)+1;progress(4)/(progress(1)+1);progress(5)/(progress(1)+1);progress(4)+acceptanceRate;progress(5)+tolScore];
%         save progress.txt -ascii progress
%         save(['temp/particle_',num2str(progress(1)),'.mat'],'models','strucIdx','xopt','tolScore','simulatedFlux');
%     end
% end
