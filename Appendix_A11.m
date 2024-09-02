%-------------------------------------------------------------------------%
init;
show_current_script_name(mfilename('fullpath'));

%-------------------------------------------------------------------------%

%% A11. Reaction time analysis: T-test comparison of RTs across Experiments 5-6 and task modalities (ES_d: trials which choose S option, ES_e: trials which choose E option, EE)

%-------------------------------------------------------------------------%
% parameters of the script                                                %
%-------------------------------------------------------------------------%
selected_exp = [5, 6];%, 6.2, 7.1, 7.2];
displayfig = 'off';
colors = [orange; green];
zscored = 0;

displayfig = 'on'
stats_data = table();

filename1 = 'Fig5B';
filename2 = 'Fig5C';

figfolder = 'fig';

figname1 = sprintf('%s/%s.svg', figfolder, filename1);
figname2 = sprintf('%s/%s.svg', figfolder, filename2);
stats_filename = sprintf('data/stats/%s.csv', [filename1 'C']);

num = 0;

lotp = [.1, .2, .3, .4, .6, .7, .8, .9];
%lotp = [0, .1, .2, .3, .4, .5,  .6, .7, .8, .9, 1];

sub_count = 0;
for exp_num = selected_exp

    num = num + 1;

    %---------------------------------------------------------------------%
    % get data parameters                                                 %
    % --------------------------------------------------------------------%
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    nsub = de.get_nsub_from_exp_num(exp_num);

    data_ed = de.extract_ES(exp_num);
    data_ee = de.extract_EE(exp_num);

    for sub = 1:nsub
        mask_lot = (ismember(data_ed.p2(sub,:), lotp));
        mask_cho1 = (data_ed.cho(sub,:)==1);
        mask_cho2 = (data_ed.cho(sub,:)==2);
        e(sub+sub_count,1) = median(data_ed.rtime(sub, logical(mask_lot.*mask_cho1)));
        d(sub+sub_count,1) = median(data_ed.rtime(sub, logical(mask_lot.*mask_cho2)));

        ee(sub+sub_count,1) = median(data_ee.rtime(sub,:));

        modalities = {'e', 'd', 'EE'};
        dd = {e(sub+sub_count,1); d(sub+sub_count,1); ee(sub+sub_count,1)};

        for mod_num = 1:3
            T1 = table(...
                sub+sub_count, exp_num, dd{mod_num},...
                {modalities{mod_num}}, 'variablenames',...
                {'subject', 'exp_num', 'RT', 'modality'}...
                );
            stats_data = [stats_data; T1];
        end
    end

    sub_count = sub_count + sub;

end

%-------------------------------------------------------------------------

raw_rt      = stats_data.RT;
raw_exp_num = stats_data.exp_num;
raw_subject = stats_data.subject;
raw_Treatment = [];
for ii=1:length(raw_rt)
    switch stats_data.modality{ii}
        case 'EE'
            raw_Treatment = [raw_Treatment; 0];
        case 'e'
            raw_Treatment = [raw_Treatment; 1];
        case 'd'
            raw_Treatment = [raw_Treatment; 2];
    end
end

stat_A11_1 = [];
%%% 1) modality_num
label1 = [0 0 2];
label2 = [2 1 1];
for ii=1:length(label1)
    y    = raw_rt;
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
    stat_A11_1 = [stat_A11_1; 0 {'modality'} label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
end
%%% 2) exp_num
label1 = [5];
label2 = [6];
for ii=1:length(label1)
    y    = raw_rt;
    x    = raw_exp_num;
    p_dd1  = y(x==label1(ii));
    p_dd2  = y(x==label2(ii));
    ss1  = raw_subject(x==label1(ii));
    ss2  = raw_subject(x==label2(ii));
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

    stat_A11_1 = [stat_A11_1; 2 {'exp_num'} label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
end

%%% 3) exp_num * modality_num
for mod_num = [0,2,1]
    y = raw_rt;
    x = raw_exp_num;
    z = raw_Treatment;
    label1 = [5];
    label2 = [6];
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
        tmp_corrected_pval = p.*(length(label1)).*length([0,1]);

        tmp_Paired = {'False'};
        tmp_Parametric = {'True'};
        tmp_P_adjusted = {'Bonf'};
        stat_A11_1 = [stat_A11_1; mod_num {'exp_num * modality'} label1(ii) label2(ii) ...
            tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval ...
            tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
    end
end

%-------------------------------------------------------------------------

stat_A11_2 = [];
y = raw_rt;
x = raw_exp_num;
z = raw_Treatment;

dd1 = y(z==0);
dd2 = y(z==1);
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

stat_A11_2 = [stat_A11_2; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A11_3 = [];
y = raw_rt;
x = raw_exp_num;
z = raw_Treatment;

dd1 = y(z==0);
dd2 = y(z==2);
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

stat_A11_3 = [stat_A11_3; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------

stat_A11_4 = [];
y = raw_rt;
x = raw_exp_num;
z = raw_Treatment;

dd1 = y(z==1);
dd2 = y(z==2);
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

stat_A11_4 = [stat_A11_4; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];

%-------------------------------------------------------------------------
stat_A11_1
stat_A11_2
stat_A11_3
stat_A11_4