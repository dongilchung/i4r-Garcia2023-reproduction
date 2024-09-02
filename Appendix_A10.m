%-------------------------------------------------------------------------
init;
%-------------------------------------------------------------------------

%% A10. Slope analysis between the inferred P(win) and the actual values: T-test comparison analysis of slopes across Experiments 6-7 and task modalities (ES, EE)
%-------------------------------------------------------------------------
selected_exp = [6, 7];

stats_filename = 'data/stats/Fig4CE.csv';

displayfig = 'on';
%-------------------------------------------------------------------------

stats_data = table();

num = 0;
sub_count = 0;
for exp_num = selected_exp
    num = num + 1;

    dED = de.extract_ES(exp_num);
    dEE = de.extract_EE(exp_num);

    corrED = mean(dED.corr,2)';
    corrEE = mean(dEE.corr,2)';
    % add ED exp_%num
    CCR{num, 1} = corrED;

    % add EE exp_%num
    CCR{num, 2} = corrEE;

    for sub = 1:dED.nsub
        modalities = {'ED', 'EE'};
        dd = {corrED, corrEE};
        for mod_num = 1:2
            T1 = table(...
                sub+sub_count, exp_num, dd{mod_num}(sub),...
                {modalities{mod_num}}, 'variablenames',...
                {'subject', 'exp_num', 'score', 'modality'}...
                );
            stats_data = [stats_data; T1];
        end
    end

    dcorr = [];
    p_sym = unique(dED.p1)';
    p_lot = unique(dED.p2)';
    for i = 1:length(p_sym)
        for j = 1:length(p_lot)
            dcorr(i, j) = ...
                (p_lot(j) < .5) * (p_sym(i) >= p_lot(j)) ...
                + (p_lot(j) > .5) * (p_sym(i) <= p_lot(j)) ;
        end
    end
    m(num) = mean(dcorr, 'all');

    sub_count = sub_count + sub;
end

%-------------------------------------------------------------------------

raw_score   = stats_data.score;
raw_exp_num = stats_data.exp_num;
raw_subject = stats_data.subject;
raw_Treatment = [];
for ii=1:length(raw_score)
    switch stats_data.modality{ii}
        case 'ED'
            raw_Treatment = [raw_Treatment; 0];
        case 'EE'
            raw_Treatment = [raw_Treatment; 1];
    end
end

stat_A10_1 = [];
%%% 1) exp_num
y    = raw_score;
x    = raw_exp_num;
p_dd1  = y(x==6);
p_dd2  = y(x==7);
ss1  = raw_subject(x==6);
ss2  = raw_subject(x==7);
dd1  = [];
for si=unique(ss1)'
    dd1 = [dd1; mean(p_dd1(si==ss1))];
end
dd2  = [];
for si=unique(ss2)'
    dd2 = [dd2; mean(p_dd2(si==ss2))];
end
[h,p,confidence_iterval,stats] = ttest2(dd1,dd2,'Alpha',0.025,'Tail','both','Vartype','unequal');

tmp_exp_num = [6 7];
tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;

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
tmp_corrected_pval = p.*(length(4));

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'False'};
tmp_Parametric = {'True'};

stat_A10_1 = [stat_A10_1; 1 {'exp_num'} 6 7 tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];

%%% 2) modality_num
label1 = [0];
label2 = [1];
for ii=1:length(label1)
    y    = raw_score;
    x    = raw_Treatment;

    p_dd1  = y(x==label1(ii));
    p_dd2  = y(x==label2(ii));
    ss1  = raw_subject(x==label1(ii));
    ss2  = raw_subject(x==label2(ii));

    dd1 = y(x==label1(ii));
    dd2 = y(x==label2(ii));
    dd = y(x==label1(ii))-y(x==label2(ii));

    [h,p,confidence_iterval,stats] = ttest(dd,0,'Alpha',0.025,'Tail','both');

    tmp_exp_num = [label1(ii) label2(ii)];
    tmp_tval = stats.tstat;
    tmp_uncorrected_pval = p;
    tmp_dof = stats.df;

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

    tmp_P_adjusted = {'Bonf'};
    tmp_alternative = {'two-sided'};
    tmp_Paired = {'True'};
    tmp_Parametric = {'True'};
    stat_A10_1 = [stat_A10_1; 2 {'modality'} label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
end

%%% 3) exp_num * modality_num
for exp_num = selected_exp
    y = raw_score;
    x = raw_exp_num;
    z = raw_Treatment;
    label1 = [1];
    label2 = [0];
    for ii=1:length(label1)
        dd1 = y(x==exp_num&z==label1(ii));
        dd2 = y(x==exp_num&z==label2(ii));
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
        tmp_corrected_pval = p.*(length(label1)).*length(selected_exp);

        tmp_Paired = {'True'};
        tmp_Parametric = {'True'};
        tmp_P_adjusted = {'Bonf'};
        stat_A10_1 = [stat_A10_1; tmp_exp_num {'exp_num * modality'} ...
            label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
    end
end

%-------------------------------------------------------------------------

stat_A10_2 = [];

%%% 1) modality_num
label1 = [0];
label2 = [1];
for ii=1:length(label1)
    y    = raw_score;
    x    = raw_Treatment;

    p_dd1  = y(x==label1(ii));
    p_dd2  = y(x==label2(ii));
    ss1  = raw_subject(x==label1(ii));
    ss2  = raw_subject(x==label2(ii));

    dd1 = y(x==label1(ii));
    dd2 = y(x==label2(ii));
    dd = y(x==label1(ii))-y(x==label2(ii));

    [h,p,confidence_iterval,stats] = ttest(dd,0,'Alpha',0.025,'Tail','both');

    tmp_exp_num = [label1(ii) label2(ii)];
    tmp_tval = stats.tstat;
    tmp_uncorrected_pval = p;
    tmp_dof = stats.df;

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

    tmp_P_adjusted = {'Bonf'};
    tmp_alternative = {'two-sided'};
    tmp_Paired = {'True'};
    tmp_Parametric = {'True'};
    stat_A10_2 = [stat_A10_2; 2 {'modality'} label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
end

%%% 2) exp_num
y    = raw_score;
x    = raw_exp_num;
p_dd1  = y(x==6);
p_dd2  = y(x==7);
ss1  = raw_subject(x==6);
ss2  = raw_subject(x==7);
dd1  = [];
for si=unique(ss1)'
    dd1 = [dd1; mean(p_dd1(si==ss1))];
end
dd2  = [];
for si=unique(ss2)'
    dd2 = [dd2; mean(p_dd2(si==ss2))];
end
[h,p,confidence_iterval,stats] = ttest2(dd1,dd2,'Alpha',0.025,'Tail','both','Vartype','unequal');

tmp_exp_num = [6 7];
tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;

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
tmp_corrected_pval = p.*(length(4));

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'False'};
tmp_Parametric = {'True'};

stat_A10_2 = [stat_A10_2; 1 {'exp_num'} 6 7 tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];

%%% 3) exp_num * modality_num
for mod_num = [0,1]
    y = raw_score;
    x = raw_exp_num;
    z = raw_Treatment;
    label1 = [6];
    label2 = [7];
    for ii=1:length(label1)
        p_dd1  = y(z==mod_num&x==label1(ii));
        p_dd2  = y(z==mod_num&x==label2(ii));
        ss1  = raw_subject(z==mod_num&x==label1(ii));
        ss2  = raw_subject(z==mod_num&x==label2(ii));
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
        tmp_cond = mod_num;
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
        tmp_corrected_pval = p.*(length([0,1]));

        tmp_Paired = {'True'};
        tmp_Parametric = {'True'};
        tmp_P_adjusted = {'Bonf'};
        stat_A10_2 = [stat_A10_2; tmp_exp_num {'exp_num * modality'} ...
            label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
    end
end

%-------------------------------------------------------------------------

stat_A10_3 = [];
y = raw_score;
x = raw_exp_num;
z = raw_Treatment;

dd1 = y(x==7&z==1);
dd2 = y(x==7&z==0);
dd = dd1-dd2;

[h,p,confidence_iterval,stats] = ttest(dd,0,'Alpha',0.025,'Tail','both');

tmp_CI98 = confidence_iterval;
tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;

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
tmp_corrected_pval = p;

%%% Power
tmp_power = sampsizepwr('t', [mean(dd) std(dd)], 0, [], N);

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'True'};
tmp_Parametric = {'True'};

stat_A10_3 = [stat_A10_3; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A10_4 = [];
y = raw_score;
x = raw_exp_num;
z = raw_Treatment;

p_dd1  = y(x==6&z==1);
p_dd2  = y(x==7&z==0);
ss1  = raw_subject(x==6&z==1);
ss2  = raw_subject(x==7&z==0);
dd1  = [];
for si=unique(ss1)'
    dd1 = [dd1; mean(p_dd1(si==ss1))];
end
dd2  = [];
for si=unique(ss2)'
    dd2 = [dd2; mean(p_dd2(si==ss2))];
end
[h,p,confidence_iterval,stats] = ttest2(dd1,dd2,'Alpha',0.025,'Tail','both','Vartype','unequal');

tmp_CI98 = confidence_iterval;
tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;

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
tmp_corrected_pval = p.*(length([1]));

%%% Power
tmp_power = sampsizepwr('t2', [mean(dd1) mean(dd2)], 0, [], N);

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'True'};
tmp_Parametric = {'True'};

stat_A10_4 = [stat_A10_4; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A10_5 = [];
y = raw_score;
x = raw_exp_num;
z = raw_Treatment;

p_dd1  = y(x==6&z==1);
p_dd2  = y(x==7&z==1);
ss1  = raw_subject(x==6&z==1);
ss2  = raw_subject(x==7&z==1);
dd1  = [];
for si=unique(ss1)'
    dd1 = [dd1; mean(p_dd1(si==ss1))];
end
dd2  = [];
for si=unique(ss2)'
    dd2 = [dd2; mean(p_dd2(si==ss2))];
end
[h,p,confidence_iterval,stats] = ttest2(dd1,dd2,'Alpha',0.025,'Tail','both','Vartype','unequal');

tmp_CI98 = confidence_iterval;
tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;

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
tmp_corrected_pval = p.*(length([1]));

%%% Power
tmp_power = sampsizepwr('t2', [mean(dd1) mean(dd2)], 0, [], N);

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'True'};
tmp_Parametric = {'True'};

stat_A10_5 = [stat_A10_5; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A10_6 = [];
y = raw_score;
x = raw_exp_num;
z = raw_Treatment;

p_dd1  = y(x==6&z==0);
p_dd2  = y(x==7&z==0);
ss1  = raw_subject(x==6&z==0);
ss2  = raw_subject(x==7&z==0);
dd1  = [];
for si=unique(ss1)'
    dd1 = [dd1; mean(p_dd1(si==ss1))];
end
dd2  = [];
for si=unique(ss2)'
    dd2 = [dd2; mean(p_dd2(si==ss2))];
end
[h,p,confidence_iterval,stats] = ttest2(dd1,dd2,'Alpha',0.025,'Tail','both','Vartype','unequal');

tmp_CI98 = confidence_iterval;
tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;

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
tmp_corrected_pval = p.*(length([1]));

%%% Power
tmp_power = sampsizepwr('t2', [mean(dd1) mean(dd2)], 0, [], N);

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'True'};
tmp_Parametric = {'True'};

stat_A10_6 = [stat_A10_6; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A10_7 = [];
y = raw_score;
x = raw_exp_num;
z = raw_Treatment;

dd1 = y(x==6&z==1);
dd2 = y(x==6&z==0);
dd = dd1-dd2;
[h,p,confidence_iterval,stats] = ttest(dd,0,'Alpha',0.025,'Tail','both');

tmp_CI98 = confidence_iterval;
tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;

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
tmp_corrected_pval = p.*(length([1]));

%%% Power
tmp_power = sampsizepwr('t', [mean(dd) std(dd)], 0, [], N);

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'True'};
tmp_Parametric = {'True'};

stat_A10_7 = [stat_A10_7; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A10_8 = [];
y = raw_score;
x = raw_exp_num;
z = raw_Treatment;

p_dd1  = y(x==6&z==0);
p_dd2  = y(x==7&z==1);
ss1  = raw_subject(x==6&z==0);
ss2  = raw_subject(x==7&z==1);
dd1  = [];
for si=unique(ss1)'
    dd1 = [dd1; mean(p_dd1(si==ss1))];
end
dd2  = [];
for si=unique(ss2)'
    dd2 = [dd2; mean(p_dd2(si==ss2))];
end
[h,p,confidence_iterval,stats] = ttest2(dd1,dd2,'Alpha',0.025,'Tail','both','Vartype','unequal');

tmp_CI98 = confidence_iterval;
tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;

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
tmp_corrected_pval = p.*(length([1]));

%%% Power
tmp_power = sampsizepwr('t2', [mean(dd1) mean(dd2)], 0, [], N);

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'True'};
tmp_Parametric = {'True'};

stat_A10_8 = [stat_A10_8; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A10_9 = [];
y = raw_score;
x = raw_exp_num;
z = raw_Treatment;

p_dd1  = y(x==7&z==1);
p_dd2  = y(x==6&z==1);
ss1  = raw_subject(x==7&z==1);
ss2  = raw_subject(x==6&z==1);
dd1  = [];
for si=unique(ss1)'
    dd1 = [dd1; mean(p_dd1(si==ss1))];
end
dd2  = [];
for si=unique(ss2)'
    dd2 = [dd2; mean(p_dd2(si==ss2))];
end
[h,p,confidence_iterval,stats] = ttest2(dd1,dd2,'Alpha',0.025,'Tail','both','Vartype','unequal');

tmp_CI98 = confidence_iterval;
tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;

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
tmp_corrected_pval = p.*(length([1]));

%%% Power
tmp_power = sampsizepwr('t2', [mean(dd1) mean(dd2)], 0, [], N);

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'True'};
tmp_Parametric = {'True'};

stat_A10_9 = [stat_A10_9; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A10_10 = [];
y = raw_score;
x = raw_exp_num;
z = raw_Treatment;

p_dd1  = y(x==7&z==0);
p_dd2  = y(x==6&z==0);
ss1  = raw_subject(x==7&z==0);
ss2  = raw_subject(x==6&z==0);
dd1  = [];
for si=unique(ss1)'
    dd1 = [dd1; mean(p_dd1(si==ss1))];
end
dd2  = [];
for si=unique(ss2)'
    dd2 = [dd2; mean(p_dd2(si==ss2))];
end
[h,p,confidence_iterval,stats] = ttest2(dd1,dd2,'Alpha',0.025,'Tail','both','Vartype','unequal');

tmp_CI98 = confidence_iterval;
tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;

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
tmp_corrected_pval = p.*(length([1]));

%%% Power
tmp_power = sampsizepwr('t2', [mean(dd1) mean(dd2)], 0, [], N);

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'True'};
tmp_Parametric = {'True'};

stat_A10_10 = [stat_A10_10; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------
stat_A10_1
stat_A10_2
stat_A10_3
stat_A10_4
stat_A10_5
stat_A10_6
stat_A10_7
stat_A10_8
stat_A10_9
stat_A10_10