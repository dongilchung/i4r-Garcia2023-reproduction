% -------------------------------------------------------------------%
% This script finds the best fitting Values for each exp             %
% then plots the option value                                        %
% -------------------------------------------------------------------%
init;
% -------------------------------------------------------------------%

selected_exp = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
%selected_exp = selected_exp(1);
sessions = [0, 1];

fit_folder = 'data/repl_fit/';

nfpm = [2, 4];

repl_alpha = [];
repl_beta  = [];
repl_exp   = [];
repl_subj  = [];
force = 1;
for exp_num = selected_exp
    exp_num
    % -------------------------------------------------------------------%
    % LEARNING
    % -------------------------------------------------------------------%
    data = de.extract_LE(exp_num);
    % set parameters
    fit_params.cho = data.cho;
    fit_params.cfcho = data.cfcho; % just reversed cho
    fit_params.out = data.out==1;
    fit_params.cfout = data.cfout==1;
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

    for subject = 1:fit_params.nsub
        subject
        global xData sData;
        xData = fit_params;
        sData = subject;
        idx = 1;
        MLE = [];
        for betaPre0 = -5:.5:5
            for alphaPre0 = -5:.5:5
                betaInit=[alphaPre0, betaPre0];
                [betaHat,fval,exitflag] = fminsearch(@RescorlaWagner_utilFunc, betaInit, optimset('MaxIter',25000,'MaxFunEvals',25000));
                alphaPost = 1/(1 + exp(-betaHat(1)));
                betaPost  = 50/(1 + exp(-betaHat(2)));

                MLE(idx,:) = [exitflag, fval, alphaPost, betaPost];

                idx=idx+1;
            end
        end
        raDatabase=[];
        temp=MLE((MLE(:,2)==min(MLE(:,2))),2:end);
        sizeTemp=size(temp,1);
        if(sizeTemp==1)
            % % second 1, original number of solution from the estimation
            raDatabase=[mean(temp,1) 1 NaN];
        else
            temp2=round(temp*10^4)/10^4; % differences smaller than 10^-4 is ignored
            for tempIndex1=1:sizeTemp-1
                for tempIndex2=tempIndex1+1:sizeTemp
                    if sum(temp2(tempIndex1,:)-temp2(tempIndex2,:))~=0 % show 10^-4 difference
                        raDatabase=[mean(temp,1) size(temp,1) NaN]; % show mulitiple local minimum
                        break
                    else
                        raDatabase=[temp2(tempIndex1,:) 1 sizeTemp];
                        % fval rho lamda
                    end
                end
            end
        end
        alpha = raDatabase(:,2);
        beta = raDatabase(:,3);
        repl_alpha = [repl_alpha; alpha];
        repl_beta  = [repl_beta; beta];
        repl_exp   = [repl_exp; exp_num];
        repl_subj  = [repl_subj; subject];
    end
end
save('data/repl_fit/mh_fit_LE','repl_*')

