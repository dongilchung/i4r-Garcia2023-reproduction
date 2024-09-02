%-------------------------------------------------------------------------
init;
%-------------------------------------------------------------------------

%% A3. Correct choice rate in the LE phase: Within-experiment regression and ANOVA analyses of performance across four experiential options (Condition 1: smallest EV difference; Condition 4: largest EV difference)

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
end

% ------------------------------------------------------------------------

stat_A3_1 = [];

for exp_num = selected_exp
    y = T_con.score;
    x = T_con.cond;
    z = T_con.exp_num;
    s = T_con.subject;

    reg_x = x(z==exp_num);
    reg_y = y(z==exp_num);
    reg_s = s(z==exp_num);
    cond1 = reg_x==1;
    cond2 = reg_x==2;
    cond3 = reg_x==3;
    cond4 = reg_x==4;

    varNames = ["score",'subject',"Cond"];
    tbl = table(reg_y,reg_s,reg_x,'VariableNames',varNames);
    formula = 'score ~ 1 + Cond';

    mdl = fitlm(tbl,formula);
    tmp_CI = coefCI(mdl,0.025);

    tmp_coef = mdl.Coefficients.Estimate(1);
    tmp_se   = mdl.Coefficients.SE(1);
    tmp_tval = mdl.Coefficients.tStat(1);
    tmp_pval = mdl.Coefficients.pValue(1);

    tmp_CI98 = tmp_CI(1,:);

    tmp_dof_model = 1;
    tmp_dof_resid = mdl.DFE;
    tmp_rsquare_ord = mdl.Rsquared.Ordinary;
    tmp_rsquare_adj = mdl.Rsquared.Adjusted;
    tmp_name = {'Intercept'};
    stat_A3_1 = [stat_A3_1; exp_num tmp_name tmp_coef tmp_se tmp_tval tmp_pval tmp_rsquare_ord tmp_rsquare_adj tmp_CI98 tmp_dof_model tmp_dof_resid];

    tmp_coef = mdl.Coefficients.Estimate(2);
    tmp_se   = mdl.Coefficients.SE(2);
    tmp_tval = mdl.Coefficients.tStat(2);
    tmp_pval = mdl.Coefficients.pValue(2);
    tmp_CI98 = tmp_CI(2,:);

    tmp_dof_model = 1;
    tmp_dof_resid = mdl.DFE;
    tmp_rsquare_ord = mdl.Rsquared.Ordinary;
    tmp_rsquare_adj = mdl.Rsquared.Adjusted;
    tmp_name = {'Cond'};
    stat_A3_1 = [stat_A3_1; exp_num tmp_name tmp_coef tmp_se tmp_tval tmp_pval tmp_rsquare_ord tmp_rsquare_adj tmp_CI98 tmp_dof_model tmp_dof_resid];
end

% ------------------------------------------------------------------------

stat_A3_2 = table;
for exp_num = selected_exp
    y = T_con.score;
    x = T_con.cond;
    z = T_con.exp_num;
    s = T_con.subject;

    rm_x = x(z==exp_num);
    rm_y = y(z==exp_num);
    rm_s = s(z==exp_num);

    subjects= rm_s(rm_x==1);
    cond1   = rm_y(rm_x==1);
    cond2   = rm_y(rm_x==2);
    cond3   = rm_y(rm_x==3);
    cond4   = rm_y(rm_x==4);

    if exp_num==4
        rm_t = table(subjects,cond1,cond4,'VariableNames',{'Subject','Cond1','Cond4'});
        rm = fitrm(rm_t, 'Cond1-Cond4 ~ 1','WithinDesign', {'Cond1','Cond4'});
    else
        rm_t = table(subjects,cond1,cond2,cond3,cond4,'VariableNames',{'Subject','Cond1','Cond2','Cond3','Cond4'});
        rm = fitrm(rm_t, 'Cond1-Cond4 ~ 1','WithinDesign', {'Cond1','Cond2','Cond3','Cond4'});
    end

    ranovatbl = ranova(rm);
    ranovatbl.Exp_num = ones(size(ranovatbl,1),1).*exp_num;

    ranovatbl.Row{1} = sprintf('(Intercept):Time_%d', exp_num);
    ranovatbl.Row{2} = sprintf('Error(Time)_%d', exp_num);

    stat_A3_2 = [stat_A3_2; ranovatbl];
end

% ------------------------------------------------------------------------

stat_A3_1
stat_A3_2