% -------------------------------------------------------------------%
% This script finds the best fitting Values for each exp             %
% then plots the option value                                        %
% -------------------------------------------------------------------%
init;
% -------------------------------------------------------------------%

selected_exp = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
%selected_exp = selected_exp(1);
sessions = [0, 1];

fit_folder = 'data/pr/';


nfpm = [2, 4];

force = 1;

for exp_num = selected_exp

    %---------------------------------------------------------------------%
    % get data parameters                                                           %
    % --------------------------------------------------------------------%
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    nsub = de.get_nsub_from_exp_num(exp_num);

    sim_params.exp_num = exp_num;
    sim_params.de = de;
    sim_params.sess = sess;
    sim_params.exp_name = name;
    sim_params.nsub = nsub;

    sim_params.model = 1;
    sim_params.path = 'data/fit/learning_LE_%s_session_%d';

    [~, throw] = get_qvalues(sim_params);

    pre.alpha1 = throw.alpha1;
    pre.beta1  = throw.beta1;
    
    fprintf('Fitting exp. %s \n', num2str(exp_num));
    
    % -------------------------------------------------------------------%
    % LEARNING
    % -------------------------------------------------------------------%
    data = de.extract_LE(exp_num);
    % set parameters
    fit_params.cho_temp = data.cho;
    fit_params.out_temp = data.out==1;
    fit_params.cfout_temp = data.cfout==1;
    fit_params.con = data.con;
    fit_params.fit_cf = (exp_num>2);
    fit_params.ntrials = size(data.cho, 2);
    fit_params.model = 1;
    fit_params.nsub = data.nsub;
    fit_params.sess = data.sess;
    fit_params.exp_num = num2str(exp_num);
    fit_params.decision_rule = 1;
    fit_params.q = 0.5;
    fit_params.noptions = 2;
    fit_params.ncond = length(unique(data.con));
    
    % simulate choice
    fit_params.cho = nan(size(data.cho,1),size(data.cho,2));
    fit_params.cfcho = nan(size(data.cho,1),size(data.cho,2));
    fit_params.out = nan(size(data.cho,1),size(data.cho,2));
    fit_params.cfout = nan(size(data.cho,1),size(data.cho,2));
    for sub = 1:fit_params.nsub
        alpha = pre.alpha1(sub);
        beta  = pre.beta1(sub);

        xCon= fit_params.con(sub,:);
        nCon= length(unique(xCon));
        roundsx   = size(xCon,2);

        xChoTemp = fit_params.cho_temp(sub,:);
        xOut= fit_params.out_temp(sub,:);
        xOutCf= fit_params.cfout_temp(sub,:); % unchosen outcome

        xOut1 = [];
        xOut2 = [];
        for trial = 1:roundsx
            if xChoTemp(trial)==1
                xOut1 = [xOut1; xOut(trial)];
                xOut2 = [xOut2; xOutCf(trial)];
            else
                xOut1 = [xOut1; xOutCf(trial)];
                xOut2 = [xOut2; xOut(trial)];
            end
        end

        Qset1  = zeros(nCon, roundsx);
        Qset2  = zeros(nCon, roundsx);
        Qset1(:,1) = 0.5;
        Qset2(:,1) = 0.5;
        for trial = 1:roundsx-1
            Qset1(:,trial+1) = Qset1(:,trial);
            Qset2(:,trial+1) = Qset2(:,trial);
            tempQ1 = Qset1(xCon(trial),trial);
            tempQ2 = Qset2(xCon(trial),trial);
            phat = 1./(1+exp(-beta*(tempQ1-tempQ2)));
            if phat>rand
                tempCh = 1;
                tempCfCh = 2;
                chosenR  = xOut1(trial);
                unchosenR= xOut2(trial);
            else
                tempCh = 2;
                tempCfCh = 1;
                chosenR  = xOut2(trial);
                unchosenR= xOut1(trial);
            end
            fit_params.cho(sub,trial) = tempCh;
            fit_params.cfcho(sub,trial) = tempCfCh;
            fit_params.out(sub,trial) = chosenR;
            fit_params.cfout(sub,trial) = unchosenR;
            if fit_params.fit_cf
                if tempCh==1
                    tempPEc = chosenR - tempQ1;
                    tempPEu = unchosenR - tempQ2;
                    tempQ1 = tempQ1 + alpha.*tempPEc;
                    tempQ2 = tempQ2 + alpha.*tempPEu;
                else
                    tempPEc = chosenR - tempQ2;
                    tempPEu = unchosenR - tempQ1;
                    tempQ1 = tempQ1 + alpha.*tempPEu;
                    tempQ2 = tempQ2 + alpha.*tempPEc;
                end
            else
                if tempCh==1
                    tempPEc = chosenR - tempQ1;
                    tempPEu = 0;
                    tempQ1 = tempQ1 + alpha.*tempPEc;
                    tempQ2 = tempQ2 + alpha.*tempPEu;
                else
                    tempPEc = chosenR - tempQ2;
                    tempPEu = 0;
                    tempQ1 = tempQ1 + alpha.*tempPEu;
                    tempQ2 = tempQ2 + alpha.*tempPEc;
                end
            end
            Qset1(xCon(trial),trial+1) = tempQ1;
            Qset2(xCon(trial),trial+1) = tempQ2;
        end
        trial = roundsx;
        tempQ1 = Qset1(xCon(trial),trial);
        tempQ2 = Qset2(xCon(trial),trial);
        phat = 1./(1+exp(-beta*(tempQ1-tempQ2)));
        if phat>rand
            tempCh = 1;
            chosenR  = xOut1(trial);
            unchosenR= xOut2(trial);
        else
            tempCh = 2;
            chosenR  = xOut2(trial);
            unchosenR= xOut1(trial);
        end
        fit_params.cho(sub,trial) = tempCh;
        fit_params.cfcho(sub,trial) = tempCfCh;
        fit_params.out(sub,trial) = chosenR;
        fit_params.cfout(sub,trial) = unchosenR;
    end

    
    save_params.fit_file = sprintf(...
        '%s%s%s%s%d', fit_folder, 'learning_LE',  data.name,  '_session_', data.sess);
    
    % fmincon params
    fmincon_params.init_value = {[1, .5], [0, .5, .5],[0, .5]};
    fmincon_params.lb = {[0.001, 0.001], [0, 0, 0], [0, 0]};
    fmincon_params.ub = {[100, 1], [100, 1, 1], [100, 1]};
    
    try
        data = load(save_params.fit_file);
        
        fit_params.params = data.data('parameters');  %% Optimization parameters
        ll = data.data('ll');
        
        if force
            error('Force = True');
        end
    catch
        [fit_params.params, ll] = runfit_learning(...
            fit_params, save_params, fmincon_params);
        
    end
    
end

    
function [parameters,ll] = ...
    runfit_learning(fit_params, save_params, fmincon_params)

   
    options = optimset(...
        'Algorithm',...
        'interior-point',...
        'Display', 'off',...
        'MaxIter', 10000,...
        'MaxFunEval', 10000);

    w = waitbar(0, 'Fitting subject');
    
    tStart = tic;
    for sub = 1:fit_params.nsub
        
        waitbar(...
            sub/fit_params.nsub,...  % Compute progression
            w,...
            sprintf('%s%d%s%s', 'Fitting subject ', sub, ' in Exp. ', fit_params.exp_num)...
            );
        
        for model = fit_params.model
         
            
            [
                p1,...
                l1,...
                rep1,...
                grad1,...
                hess1,...
            ] = fmincon(...
                @(x) getlpp_learning(...
                    x,...
                    fit_params.con(sub, :),...
                    fit_params.cho(sub, :),...
                    fit_params.cfcho(sub, :),...
                    fit_params.out(sub, :),...
                    fit_params.cfout(sub, :),...
                    fit_params.q,...
                    fit_params.ntrials, model, fit_params.decision_rule,...
                    fit_params.fit_cf),...
                fmincon_params.init_value{model},...
                [], [], [], [],...
                fmincon_params.lb{model},...
                fmincon_params.ub{model},...
                [],...
                options...
                );
            
            parameters{model}(sub, :) = p1;
            ll(model, sub) = l1;

        end
    end
   toc(tStart);
    % Save the data
   %data = load(save_params.fit_file);
      
   %hessian = data.data('hessian');
   data = containers.Map({'parameters', 'll'},...
            {parameters, ll});
   save(save_params.fit_file, 'data');
     close(w);
%     
end
