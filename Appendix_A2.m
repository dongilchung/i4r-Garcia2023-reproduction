%-------------------------------------------------------------------------
init;
%-------------------------------------------------------------------------

%% A2. Correct choice rate comparison in the LE phase: Pairwise t-tests of performance across Experiments 1-3 and four experiential options (Condition 1: smallest EV difference; Condition 4: largest EV difference)

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

stat_A2_1 = [];
% All
label1 = [1 1 2];
label2 = [2 3 3];
for ii=1:length(label1)
    y    = T_con.score;
    x    = T_con.exp_num;
    p_dd1  = y(x==label1(ii));
    p_dd2  = y(x==label2(ii));
    ss1  = T_con.subject(x==label1(ii));
    ss2  = T_con.subject(x==label2(ii));
    dd1  = [];
    for si=unique(ss1)'
        dd1 = [dd1; mean(p_dd1(si==ss1))];
    end
    dd2  = [];
    for si=unique(ss2)'
        dd2 = [dd2; mean(p_dd2(si==ss2))];
    end

    [h,p,confidence_iterval,stats] = ttest2(dd1,dd2,'Alpha',0.025,'Tail','both','Vartype','unequal');

    tmp_exp_num = [label1(ii) label2(ii)];
    tmp_cond = 0;
    tmp_tval = stats.tstat;
    tmp_uncorrected_pval = p;
    tmp_dof = stats.df;
    tmp_CI98 = (confidence_iterval);
    tmp_alternative = {'two-sided'};

    %%% Bayes Factor
    T = stats.tstat;
    df = stats.df;
    N = prod([numel(dd1) numel(dd2)])/sum([numel(dd1) numel(dd2)]);
    r = sqrt(2)/2;
    numerator = (1+T.^2/df).^(-(df+1)/2);
    fun  = @(g) ( ((1+N.*g.*r.^2).^-0.5) .* (1+T.^2./((1+N.*g.*r.^2).*df)).^(-(df+1)/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );

    bf01 = numerator/integral(fun,0,inf);
    bf10 = 1./bf01;
    tmp_BF10 = bf10;

    %%% Cohen-d
    tmp_cohen_d = (mean(dd1)-mean(dd2))./sqrt((std(dd1).^2+std(dd2).^2)/2);

    %%% corrected-p
    tmp_corrected_pval = p.*(length(label1));

    tmp_Paired = {'False'};
    tmp_Parametric = {'True'};
    tmp_P_adjusted = {'Bonf'};

    stat_A2_1 = [stat_A2_1; 0 label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
end

% ------------------------------------------------------------------------

label1 = [1 1 1 2 2 3];
label2 = [2 3 4 3 4 4];
for ii=1:length(label1)
    y    = T_con.score;
    x    = T_con.cond;
    z    = T_con.exp_num;
    x(z==4)=[];
    y(z==4)=[];

    dd1 = y(x==label1(ii));
    dd2 = y(x==label2(ii));
    dd = y(x==label1(ii))-y(x==label2(ii));

    [h,p,confidence_iterval,stats] = ttest(dd,0,'Alpha',0.025,'Tail','both');

    tmp_exp_num = 0;
    tmp_cond = [label1(ii) label2(ii)];
    tmp_tval = stats.tstat;
    tmp_uncorrected_pval = p;
    tmp_dof = stats.df;
    tmp_CI98 = (confidence_iterval);
    tmp_alternative = {'two-sided'};

    %%% Bayes Factor
    T = stats.tstat;
    df = stats.df;
    N = numel(dd);
    r = sqrt(2)/2;
    numerator = (1+T.^2/df).^(-(df+1)/2);
    fun  = @(g) ( ((1+N.*g.*r.^2).^-0.5) .* (1+T.^2./((1+N.*g.*r.^2).*df)).^(-(df+1)/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );

    bf01 = numerator/integral(fun,0,inf);
    bf10 = 1./bf01;
    tmp_BF10 = bf10;

    %%% Cohen-d
    tmp_cohen_d = (mean(dd1)-mean(dd2))./sqrt((std(dd1).^2+std(dd2).^2)/2);

    %%% corrected-p
    tmp_corrected_pval = p.*(length(label1));

    tmp_Paired = {'True'};
    tmp_Parametric = {'True'};
    tmp_P_adjusted = {'Bonf'};

    stat_A2_1 = [stat_A2_1; 0 label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
end

%-------------------------------------------------------------------------

for exp_num = selected_exp(1:3)
    y = T_con.score;
    x = T_con.cond;
    z = T_con.exp_num;

    label1 = [1 1 1 2 2 3];
    label2 = [2 3 4 3 4 4];
    for ii=1:length(label1)
        dd1 = y(z==exp_num&x==label1(ii));
        dd2 = y(z==exp_num&x==label2(ii));
        dd  = dd1-dd2;
        [h,p,confidence_iterval,stats] = ttest(dd,0,'Alpha',0.025,'Tail','both');
        tmp_exp_num = exp_num;
        tmp_cond = [label1(ii) label2(ii)];
        tmp_tval = stats.tstat;
        tmp_uncorrected_pval = p;
        tmp_dof = stats.df;
        tmp_CI98 = (confidence_iterval);
        tmp_alternative = {'two-sided'};

        %%% Bayes Factor
        T = stats.tstat;
        df = stats.df;
        N = numel(dd);
        r = sqrt(2)/2;
        numerator = (1+T.^2/df).^(-(df+1)/2);
        fun  = @(g) ( ((1+N.*g.*r.^2).^-0.5) .* (1+T.^2./((1+N.*g.*r.^2).*df)).^(-(df+1)/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );

        bf01 = numerator/integral(fun,0,inf);
        bf10 = 1./bf01;
        tmp_BF10 = bf10;

        %%% Cohen-d
        tmp_cohen_d = (mean(dd1)-mean(dd2))./sqrt((std(dd1).^2+std(dd2).^2)/2);

        %%% corrected-p
        tmp_corrected_pval = p.*(length(label1)).*length(selected_exp(1:3));

        tmp_Paired = {'True'};
        tmp_Parametric = {'True'};
        tmp_P_adjusted = {'Bonf'};
        stat_A2_1 = [stat_A2_1; tmp_exp_num label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
    end
end

%-------------------------------------------------------------------------

stat_A2_2 = [];
y = T_con.score;
x = T_con.cond;
z = T_con.exp_num;

dd1 = y(z==1&x==1);
dd2 = y(z==1&x==4);
dd  = dd1-dd2;
[h,p,confidence_iterval,stats] = ttest(dd,0,'Alpha',0.025,'Tail','both');
tmp_exp_num = 1;
tmp_tval = stats.tstat;
tmp_pval = p;
tmp_dof = stats.df;
tmp_CI98 = (confidence_iterval);
tmp_alternative = {'two-sided'};

%%% Bayes Factor
T = stats.tstat;
df = stats.df;
N = numel(dd);
r = sqrt(2)/2;
numerator = (1+T.^2/df).^(-(df+1)/2);
fun  = @(g) ( ((1+N.*g.*r.^2).^-0.5) .* (1+T.^2./((1+N.*g.*r.^2).*df)).^(-(df+1)/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );

bf01 = numerator/integral(fun,0,inf);
bf10 = 1./bf01;
tmp_BF10 = bf10;

%%% Cohen-d
tmp_cohen_d = (mean(dd1)-mean(dd2))./sqrt((std(dd1).^2+std(dd2).^2)/2);

%%% Power
tmp_power = sampsizepwr('t', [mean(dd) std(dd)], 0, [], N);

stat_A2_2 = [stat_A2_2; tmp_tval tmp_dof tmp_alternative tmp_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A2_3 = [];
y = T_con.score;
x = T_con.cond;
z = T_con.exp_num;

dd1 = y(z==2&x==1);
dd2 = y(z==2&x==4);
dd  = dd1-dd2;
[h,p,confidence_iterval,stats] = ttest(dd,0,'Alpha',0.025,'Tail','both');
tmp_exp_num = 1;
tmp_tval = stats.tstat;
tmp_pval = p;
tmp_dof = stats.df;
tmp_CI98 = (confidence_iterval);
tmp_alternative = {'two-sided'};

%%% Bayes Factor
T = stats.tstat;
df = stats.df;
N = numel(dd);
r = sqrt(2)/2;
numerator = (1+T.^2/df).^(-(df+1)/2);
fun  = @(g) ( ((1+N.*g.*r.^2).^-0.5) .* (1+T.^2./((1+N.*g.*r.^2).*df)).^(-(df+1)/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );

bf01 = numerator/integral(fun,0,inf);
bf10 = 1./bf01;
tmp_BF10 = bf10;

%%% Cohen-d
tmp_cohen_d = (mean(dd1)-mean(dd2))./sqrt((std(dd1).^2+std(dd2).^2)/2);

%%% Power
tmp_power = sampsizepwr('t', [mean(dd) std(dd)], 0, [], N);

stat_A2_3 = [stat_A2_3; tmp_tval tmp_dof tmp_alternative tmp_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A2_1
stat_A2_2
stat_A2_3