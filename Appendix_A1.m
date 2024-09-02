%-------------------------------------------------------------------------
init;
%-------------------------------------------------------------------------

%% A1. Correct choice rate analysis in the LE phase: Comparing individual performance against chance level (0.5) using t-tests across Experiments 1-4

% ------------------------------------------------------------------------
% customizable parameters
% ------------------------------------------------------------------------
selected_exp = [1, 2, 3, 4];
modality = 'LE';
color = blue;
displayfig = 'on';
filename = 'Fig2A';

% ------------------------------------------------------------------------
% fixed parameters
% ------------------------------------------------------------------------
T_exp = table();
T_con = table();


sub_count = 0;
num = 0;
stat_A1 = [];
for exp_num = selected_exp
    clear dd
    num = num + 1;

    data = de.extract_LE(exp_num);
    data_ed = de.extract_ES(exp_num);

    if exp_num == 4
        data.con(data.con == 2) = 4;
    end
    ncon = length(unique(data.con));

    dd = NaN(ncon, data.nsub);
    cons = flip(unique(data.con));

    for i = 1:ncon
        for sub = 1:data.nsub

            dd(i, sub) = mean(...
                data.corr(sub, data.con(sub,:)==cons(i)));
            %if ismember(cons(i), [1, 4])
            complete = ismember(exp_num, [3, 4]);
            block = ismember(exp_num, [2, 3, 4]);
            less_cues = exp_num == 4;


            T3 = table(...
                sub+sub_count, exp_num,  complete,  block, less_cues, dd(i, sub), cons(i), ...
                'variablenames',...
                {'subject', 'exp_num', 'complete', 'block','less_cues', 'score', 'cond'}...
                );

            T_con = [T_con; T3];
            % end
        end
    end

    for sub = 1:data.nsub
        s1 = mean(data.corr(sub, :));
        s2 = mean(data_ed.corr(sub, :));

        complete = int32(ismember(exp_num, [3, 4]));
        block = int32(ismember(exp_num, [2, 3, 4]));
        less_cues = int32(exp_num == 4);

        T1 = table(...
            sub+sub_count, exp_num, complete,  block, less_cues, s1, {'LE'}, ...
            'variablenames',...
            {'subject', 'exp_num', 'complete', 'block','less_cues', 'score', 'modality'}...
            );

        T_exp = [T_exp; T1];
    end
    addpath(genpath('./repl_utils'))
    idx = 1;
    if exp_num~=4
        num_cond = [1:4];
    else
        num_cond = [1,4];
    end
    for ci = size(dd,1):-1:1
        data = dd(ci,:);
        chance_lv = 0.5;
        [h,p,confidence_iterval,stats] = ttest(data,chance_lv,'Alpha',0.025,'Tail','both');
        % p1out = sampsizepwr('t',[stats.tstat stats.sd],[],0.98,length(data),'Tail','both')

        tmp_exp_num = exp_num;
        tmp_cond = num_cond(idx);
        tmp_tval = stats.tstat;
        tmp_pval = p;
        tmp_dof = stats.df;
        tmp_CI98 = (confidence_iterval);
        tmp_alternative = {'two-sided'};

        %%% Bayes Factor
        T = stats.tstat;
        df = stats.df;
        N = numel(data);
        r = sqrt(2)/2;
        numerator = (1+T.^2/df).^(-(df+1)/2);
        fun  = @(g) ( ((1+N.*g.*r.^2).^-0.5) .* (1+T.^2./((1+N.*g.*r.^2).*df)).^(-(df+1)/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );

        bf01 = numerator/integral(fun,0,inf);
        bf10 = 1./bf01;
        tmp_BF10 = bf10;

        %%% Cohen-d
        tmp_cohen_d = (mean(data)-0.5)./std(data);

        %%% Power
        tmp_power = sampsizepwr('t', [mean(data) std(data)], 0.5, [], N);

        stat_A1 = [stat_A1; tmp_exp_num tmp_cond tmp_tval tmp_dof tmp_alternative tmp_pval {tmp_CI98} tmp_cohen_d tmp_BF10 tmp_power];
        idx = idx+1;
    end
end

stat_A1