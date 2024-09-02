%-------------------------------------------------------------------------
init;
%-------------------------------------------------------------------------

%% A4. Correct choice rate in the LE phase: Regression analysis of performance across Experiments 1-4 and four experiential options (Condition 1: smallest EV difference; Condition 4: largest EV difference)

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

stat_A4_1 = []; % mixed LM
exp_num = 1
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

varNames = ["score",'subject',"T1","T2","T3","T4"];
tbl = table(reg_y,reg_s,cond4,cond3,cond2,cond1,'VariableNames',varNames);
formula = 'score ~ 1 + T2 + T3 + T4 + (1 + T2 + T3 + T4|subject)';

mdl = fitglme(tbl,formula);
tmp_CI = coefCI(mdl,'alpha',0.025);

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
stat_A4_1 = [stat_A4_1; exp_num tmp_name tmp_coef tmp_se tmp_tval tmp_pval tmp_rsquare_ord tmp_rsquare_adj tmp_CI98 tmp_dof_model tmp_dof_resid];

tmp_coef = mdl.Coefficients.Estimate(2);
tmp_se   = mdl.Coefficients.SE(2);
tmp_tval = mdl.Coefficients.tStat(2);
tmp_pval = mdl.Coefficients.pValue(2);
tmp_CI98 = tmp_CI(2,:);

tmp_dof_model = 1;
tmp_dof_resid = mdl.DFE;
tmp_rsquare_ord = mdl.Rsquared.Ordinary;
tmp_rsquare_adj = mdl.Rsquared.Adjusted;
tmp_name = {'C(cond)[T.2]'};
stat_A4_1 = [stat_A4_1; exp_num tmp_name tmp_coef tmp_se tmp_tval tmp_pval tmp_rsquare_ord tmp_rsquare_adj tmp_CI98 tmp_dof_model tmp_dof_resid];

tmp_coef = mdl.Coefficients.Estimate(3);
tmp_se   = mdl.Coefficients.SE(3);
tmp_tval = mdl.Coefficients.tStat(3);
tmp_pval = mdl.Coefficients.pValue(3);
tmp_CI98 = tmp_CI(3,:);

tmp_dof_model = 1;
tmp_dof_resid = mdl.DFE;
tmp_rsquare_ord = mdl.Rsquared.Ordinary;
tmp_rsquare_adj = mdl.Rsquared.Adjusted;
tmp_name = {'C(cond)[T.3]'};
stat_A4_1 = [stat_A4_1; exp_num tmp_name tmp_coef tmp_se tmp_tval tmp_pval tmp_rsquare_ord tmp_rsquare_adj tmp_CI98 tmp_dof_model tmp_dof_resid];

tmp_coef = mdl.Coefficients.Estimate(4);
tmp_se   = mdl.Coefficients.SE(4);
tmp_tval = mdl.Coefficients.tStat(4);
tmp_pval = mdl.Coefficients.pValue(4);
tmp_CI98 = tmp_CI(4,:);

tmp_dof_model = 1;
tmp_dof_resid = mdl.DFE;
tmp_rsquare_ord = mdl.Rsquared.Ordinary;
tmp_rsquare_adj = mdl.Rsquared.Adjusted;
tmp_name = {'C(cond)[T.4]'};
stat_A4_1 = [stat_A4_1; exp_num tmp_name tmp_coef tmp_se tmp_tval tmp_pval tmp_rsquare_ord tmp_rsquare_adj tmp_CI98 tmp_dof_model tmp_dof_resid];

% ------------------------------------------------------------------------

stat_A4_2 = [];

y = T_con.score;
x = T_con.cond;
z = T_con.exp_num;
s = T_con.subject;

reg_x = [];
reg_y = [];
reg_z = [];
reg_s = [];
for exp_num = selected_exp
    selected_subj = unique(s(z==exp_num));
    for sub_num = selected_subj'
        reg_x = [reg_x; mean(x(s==sub_num&z==exp_num))];
        reg_y = [reg_y; mean(y(s==sub_num&z==exp_num))];
        reg_z = [reg_z; mean(z(s==sub_num&z==exp_num))];
        reg_s = [reg_s; mean(s(s==sub_num&z==exp_num))];
    end
end

cond1 = reg_z==1;
cond2 = reg_z==2;
cond3 = reg_z==3;
cond4 = reg_z==4;

varNames = ["score",'subject',"T1","T2","T3","T4"];
tbl = table(reg_y,reg_s,cond1,cond2,cond3,cond4,'VariableNames',varNames);
formula = 'score ~ 1 + T2 + T3 + T4';

mdl = fitlm(tbl,formula,'RobustOpts','ols');
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
stat_A4_2 = [stat_A4_2; exp_num tmp_name tmp_coef tmp_se tmp_tval tmp_pval tmp_rsquare_ord tmp_rsquare_adj tmp_CI98 tmp_dof_model tmp_dof_resid];

tmp_coef = mdl.Coefficients.Estimate(2);
tmp_se   = mdl.Coefficients.SE(2);
tmp_tval = mdl.Coefficients.tStat(2);
tmp_pval = mdl.Coefficients.pValue(2);
tmp_CI98 = tmp_CI(2,:);

tmp_dof_model = 1;
tmp_dof_resid = mdl.DFE;
tmp_rsquare_ord = mdl.Rsquared.Ordinary;
tmp_rsquare_adj = mdl.Rsquared.Adjusted;
tmp_name = {'C(exp_num)[T.2]'};
stat_A4_2 = [stat_A4_2; exp_num tmp_name tmp_coef tmp_se tmp_tval tmp_pval tmp_rsquare_ord tmp_rsquare_adj tmp_CI98 tmp_dof_model tmp_dof_resid];

tmp_coef = mdl.Coefficients.Estimate(3);
tmp_se   = mdl.Coefficients.SE(3);
tmp_tval = mdl.Coefficients.tStat(3);
tmp_pval = mdl.Coefficients.pValue(3);
tmp_CI98 = tmp_CI(3,:);

tmp_dof_model = 1;
tmp_dof_resid = mdl.DFE;
tmp_rsquare_ord = mdl.Rsquared.Ordinary;
tmp_rsquare_adj = mdl.Rsquared.Adjusted;
tmp_name = {'C(exp_num)[T.3]'};
stat_A4_2 = [stat_A4_2; exp_num tmp_name tmp_coef tmp_se tmp_tval tmp_pval tmp_rsquare_ord tmp_rsquare_adj tmp_CI98 tmp_dof_model tmp_dof_resid];

tmp_coef = mdl.Coefficients.Estimate(4);
tmp_se   = mdl.Coefficients.SE(4);
tmp_tval = mdl.Coefficients.tStat(4);
tmp_pval = mdl.Coefficients.pValue(4);
tmp_CI98 = tmp_CI(4,:);

tmp_dof_model = 1;
tmp_dof_resid = mdl.DFE;
tmp_rsquare_ord = mdl.Rsquared.Ordinary;
tmp_rsquare_adj = mdl.Rsquared.Adjusted;
tmp_name = {'C(exp_num)[T.4]'};
stat_A4_2 = [stat_A4_2; exp_num tmp_name tmp_coef tmp_se tmp_tval tmp_pval tmp_rsquare_ord tmp_rsquare_adj tmp_CI98 tmp_dof_model tmp_dof_resid];

% ------------------------------------------------------------------------

stat_A4_1
stat_A4_2