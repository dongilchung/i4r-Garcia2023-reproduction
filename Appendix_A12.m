%-------------------------------------------------------------------------%
init;
show_current_script_name(mfilename('fullpath'));

%-------------------------------------------------------------------------%

%% A12. 

%-------------------------------------------------------------------------%
% parameters of the script                                                %
%-------------------------------------------------------------------------%
selected_exp = [1, 2, 3, 4, 5, 6];
displayfig = 'on';
colors = [red; dark_blue; pink; black];

%filenames
filename = 'Fig5D';
figfolder = 'fig';

figname = sprintf('%s/%s.svg', figfolder, filename);
stats_filename = sprintf('data/stats/%s.csv', filename);


% %  
% figure('Renderer', 'Painter', 'Units', 'centimeters',...
%     'Position', [0,0,5.3, 5.3/1.25], 'visible', displayfig)

sub_count = 0;
stats_data = table();
num = 0;

for exp_num = selected_exp
    
    num = num + 1;
    
    
    %---------------------------------------------------------------------%
    % get data parameters                                                 %
    % --------------------------------------------------------------------%
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    nsub = de.get_nsub_from_exp_num(exp_num);
    
    data = de.extract_ES(exp_num);
    symp = unique(data.p1);
    heur = heuristic(data);
    le = [];
    
    % get le q values estimates
    for i = 1:length(sess)
        sim_params.de = de;
        sim_params.sess = sess(i);
        sim_params.exp_name = name;
        sim_params.exp_num = exp_num;
        sim_params.nsub = nsub;
        sim_params.model = 1;
        sim_params.path = 'data/fit/learning_LE_%s_session_%d';

        
        if length(sess) == 2
            d = de.extract_ES(str2num(sprintf('%d.%d', exp_num, sess(i)+1)));
        else
            d = data;
        end
        
        [Q, tt] = get_qvalues(sim_params);

        le = [le argmax_estimate(d, symp, Q)];
        
    end
    
    for sub = 1:nsub
        s = sub + sub_count;
        
        S1 = mean(logical((...
                data.cho(sub,:)== heur(sub,:)) .* (data.cho(sub,:)~=le(sub,:))));
        S2 = mean(logical((...
                data.cho(sub,:)~= heur(sub,:)) .* (data.cho(sub,:)==le(sub,:))));
           
        o_heur(s,1) = median(...
            data.rtime(sub, logical((...
                data.cho(sub,:)== heur(sub,:)) .* (data.cho(sub,:)~=le(sub,:)))));
        o_le(s,1) = median(...
            data.rtime(sub,logical(...
            (data.cho(sub,:)~= heur(sub,:)) .* (data.cho(sub,:)==le(sub,:)))));
        
        none(s,1) = median(...
            data.rtime(sub,logical(...
            (data.cho(sub,:)~=heur(sub,:)).*(data.cho(sub,:)~=le(sub,:)))));
        both(s,1) = median(...
            data.rtime(sub,logical(...
            (data.cho(sub,:)==heur(sub,:)).*(data.cho(sub,:)==le(sub,:)))));
        
        modalities = {'heur', 'le', 'both', 'none'};
        
        %dd = {o_heur(s,1); o_le(s,1);...
        %    both(s,1); none(s,1)};
        
        %score = {S1, S2, NaN, NaN};
        rt2 = median(data.rtime(sub, :));
        
        %for mod_num = 1:4
            T1 = table(...
                s, exp_num, o_heur(s,1), {'H'}, ...
                 'variablenames',...
                {'subject', 'exp_num', 'RT', 'modality'}...
                );
            stats_data = [stats_data; T1];
      
            T1 = table(...
                s, exp_num, o_le(s,1), {'LE'}, ...
                 'variablenames',...
                {'subject', 'exp_num', 'RT', 'modality'}...
                );
             stats_data = [stats_data; T1];

    end
    
    sub_count = sub_count + sub;
    
  
end

%-------------------------------------------------------------------------%

raw_rt   = stats_data.RT;
raw_exp_num = stats_data.exp_num;
raw_subject = stats_data.subject;
raw_Treatment = [];
for ii=1:length(raw_rt)
    switch stats_data.modality{ii}
        case 'H'
            raw_Treatment = [raw_Treatment; 0];
        case 'LE'
            raw_Treatment = [raw_Treatment; 1];
    end
end

stat_A12_1 = [];
y = raw_rt;
x = raw_exp_num;
z = raw_Treatment;

dd1 = y(z==0);
dd2 = y(z==1);
dd = dd1-dd2;

[h,p,confidence_iterval,stats] = ttest(dd,0,'Alpha',0.025,'Tail','both');

tmp_tval = stats.tstat;
tmp_uncorrected_pval = p;
tmp_dof = stats.df;
tmp_CI98 = (confidence_iterval);

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
tmp_cohen_d = -(nanmean(dd1)-nanmean(dd2))./sqrt((nanstd(dd1).^2+nanstd(dd2).^2)/2);

%%% corrected-p
tmp_corrected_pval = p.*30;

%%% Power
tmp_power = sampsizepwr('t', [nanmean(dd) nanstd(dd)], 0, [], N);

tmp_P_adjusted = {'Bonf'};
tmp_alternative = {'two-sided'};
tmp_Paired = {'True'};
tmp_Parametric = {'True'};

stat_A12_1 = [stat_A12_1; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power tmp_corrected_pval];

%-------------------------------------------------------------------------%

stat_A12_2 = [];
%%% 1) modality_num
label1 = [0];
label2 = [1];
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
    tmp_cohen_d = (nanmean(dd1)-nanmean(dd2))./sqrt((nanstd(dd1).^2+nanstd(dd2).^2)/2);

    %%% corrected-p
    tmp_corrected_pval = p.*30;

    tmp_P_adjusted = {'Bonf'};
    tmp_alternative = {'two-sided'};
    tmp_Paired = {'True'};
    tmp_Parametric = {'True'};
    stat_A12_2 = [stat_A12_2; 0 {'modality'} label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
end
%%% 2) exp_num
label1 = [1 1 1 1 1 2 2 2 2 3 3 3 4 4 5];
label2 = [2 3 4 5 6 3 4 5 6 4 5 6 5 6 6];
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
    tmp_cohen_d = (nanmean(dd1)-nanmean(dd2))./sqrt((nanstd(dd1).^2+nanstd(dd2).^2)/2);

    %%% corrected-p
    tmp_corrected_pval = p.*30;

    tmp_P_adjusted = {'Bonf'};
    tmp_alternative = {'two-sided'};
    tmp_Paired = {'False'};
    tmp_Parametric = {'True'};

    stat_A12_2 = [stat_A12_2; 0 {'exp_num'} label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
end

%%% 3) exp_num * modality_num
for mod_num = [0,1]
    y = raw_rt;
    x = raw_exp_num;
    z = raw_Treatment;
    label1 = [1 1 1 1 1 2 2 2 2 3 3 3 4 4 5];
    label2 = [2 3 4 5 6 3 4 5 6 4 5 6 5 6 6];
    for ii=1:length(label1)
        p_dd1  = y(z==mod_num&x==label1(ii));
        p_dd2  = y(z==mod_num&x==label2(ii));
        ss1  = raw_subject(z==mod_num&x==label1(ii));
        ss2  = raw_subject(z==mod_num&x==label2(ii));
        dd1  = [];
        for si=unique(ss1)'
            if sum(~isnan(y(raw_subject==si)))==2
                dd1 = [dd1; mean(p_dd1(si==ss1))];
            else
                dd1 = [dd1; nan];
            end
        end
        dd2  = [];
        for si=unique(ss2)'
            if sum(~isnan(y(raw_subject==si)))==2
                dd2 = [dd2; mean(p_dd2(si==ss2))];
            else
                dd2 = [dd2; nan];
            end
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
        tmp_cohen_d = (nanmean(dd1)-nanmean(dd2))./sqrt((nanstd(dd1).^2+nanstd(dd2).^2)/2);

        %%% corrected-p
        tmp_corrected_pval = p.*30;

        tmp_Paired = {'False'};
        tmp_Parametric = {'True'};
        tmp_P_adjusted = {'Bonf'};
        stat_A12_2 = [stat_A12_2; mod_num {'exp_num * modality'} label1(ii) label2(ii) ...
            tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval ...
            tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
    end
end

%-------------------------------------------------------------------------%
stat_A12_1
stat_A12_2