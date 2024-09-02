%-------------------------------------------------------------------------
init;
show_current_script_name(mfilename('fullpath'));
%-------------------------------------------------------------------------

%% A7. Slope analysis between the inferred P(win) and the actual values: T-test comparison of slopes across Experiments 1-4 and task modalities (LE, ES, and SP)

%-------------------------------------------------------------------------%
% parameters of the script                                                %
%-------------------------------------------------------------------------%
selected_exp = [1, 2, 3, 4];
modalities = {'LE', 'ES', 'SP'};
displayfig = 'on';
colors = [blue;orange;magenta];
% filenames
filename = 'Fig2D';
figfolder = 'fig';

figname = sprintf('%s/%s.svg', figfolder, filename);
stats_filename = sprintf('data/stats/%s.csv', filename);


%-------------------------------------------------------------------------%
% prepare data                                                            %
%-------------------------------------------------------------------------%
% stats_data is table that is used to compute stats later
stats_data = table();

num = 0;
sub_count = 0;
for exp_num = selected_exp
    num = num + 1;
    disp(num)
    
    %---------------------------------------------------------------------%
    % get data parameters                                                           %
    % --------------------------------------------------------------------%
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    nsub = de.get_nsub_from_exp_num(exp_num);
    
    throw = de.extract_ES(exp_num);
    nsym = length(unique(throw.p1));
    p1 = unique(throw.p1)'.*100;
    
    % prepare data structure
    midpoints = nan(length(modalities), nsub, nsym);
    slope = nan(length(modalities), nsub, 2);
    reshape_midpoints = nan(nsub, nsym);
    
    sim_params.exp_num = exp_num;
    sim_params.de = de;
    sim_params.sess = sess;
    sim_params.exp_name = name;
    sim_params.nsub = nsub;

    
    for mod_num = 1:length(modalities)
        
        % get data depending on chosen modality
        switch (modalities{mod_num})
            
            case 'LE'
                sim_params.model = 1;
                sim_params.path = 'data/fit/learning_LE_%s_session_%d';

                [midpoints(mod_num, :, :), throw] = get_qvalues(sim_params);
                
            case {'EE', 'ES'}
               
                param = load(...
                    sprintf('data/fit/midpoints_%s_%s_session_%d',...
                    modalities{mod_num}, name, sess));

                midpoints(mod_num, :, :) = param.midpoints;
                
            case 'SP'
                sim_params.model = 2;

                [midpoints(mod_num, :, :), throw] = get_qvalues(sim_params);
        end
                  
        % fill data
        reshape_midpoints(:, :) = midpoints(mod_num, :, :);
        slope(mod_num,:,:) = add_linear_reg(...
            reshape_midpoints.*100, p1, colors(mod_num, :));
        
        % fill data for stats
        for sub = 1:nsub
            T1 = table(...
                sub+sub_count, num, slope(mod_num, sub, 2),...
                {modalities{mod_num}}, 'variablenames',...
                {'subject', 'exp_num', 'slope', 'modality'}...
                );
            stats_data = [stats_data; T1];
        end
    end
    sub_count = sub_count+sub;
end

%-------------------------------------------------------------------------

raw_slope   = stats_data.slope;
raw_exp_num = stats_data.exp_num;
raw_subject = stats_data.subject;
raw_Treatment = [];
for ii=1:length(stats_data.slope)
    switch stats_data.modality{ii}
        case 'LE'
            raw_Treatment = [raw_Treatment; 1];
        case 'ES'
            raw_Treatment = [raw_Treatment; 2];
        case 'SP'
            raw_Treatment = [raw_Treatment; 3];
    end
end

stat_A7_1 = [];
%%% 1) exp_num
label1 = [1 1 1 2 2 3];
label2 = [2 3 4 3 4 4];
for ii=1:length(label1)
    y    = raw_slope;
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
    
    tmp_exp_num = [label1(ii) label2(ii)];
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
    tmp_corrected_pval = p.*(length(label1));

    tmp_P_adjusted = {'Bonf'};
    tmp_alternative = {'two-sided'};
    tmp_Paired = {'False'};
    tmp_Parametric = {'True'};
    stat_A7_1 = [stat_A7_1; ii {'exp_num'} label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
end

%%% 2) modality_num
label1 = [2 2 1];
label2 = [1 3 3];
for ii=1:length(label1)
    y    = raw_slope;
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
    stat_A7_1 = [stat_A7_1; ii {'modality'} label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
end

%%% 3) exp_num * modality_num
for exp_num = selected_exp
    y = raw_slope;
    x = raw_exp_num;
    z = raw_Treatment;
    label1 = [2 2 1];
    label2 = [1 3 3];
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
        stat_A7_1 = [stat_A7_1; tmp_exp_num {'exp_num * modality'} ...
            label1(ii) label2(ii) tmp_Paired tmp_Parametric tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_corrected_pval tmp_P_adjusted tmp_BF10 tmp_cohen_d];
    end
end

%-------------------------------------------------------------------------

stat_A7_2 = [];
for exp_num = selected_exp
    y = raw_slope;
    x = raw_exp_num;
    z = raw_Treatment;

    p_dd1  = y(x==exp_num&z==2);
    p_dd2  = y(x==exp_num&z==3);
    ss1  = raw_subject(x==label1(ii));
    ss2  = raw_subject(x==label2(ii));
    
    dd1 = y(x==exp_num&z==2);
    dd2 = y(x==exp_num&z==3);
    dd = dd1-dd2;

    [h,p,confidence_iterval,stats] = ttest(dd,0,'Alpha',0.025,'Tail','both');
    
    tmp_exp_num = [label1(ii) label2(ii)];
    tmp_tval = stats.tstat;
    tmp_CI98 = (confidence_iterval);
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

    %%% Power
    tmp_power = sampsizepwr('t', [mean(dd) std(dd)], 0, [], N);

    tmp_P_adjusted = {'Bonf'};
    tmp_alternative = {'two-sided'};
    tmp_Paired = {'True'};
    tmp_Parametric = {'True'};
    stat_A7_2 = [stat_A7_2; {'T-test'} tmp_tval tmp_dof tmp_alternative tmp_uncorrected_pval tmp_CI98 tmp_cohen_d tmp_BF10 tmp_power];
end

%-------------------------------------------------------------------------

stat_A7_1
stat_A7_2