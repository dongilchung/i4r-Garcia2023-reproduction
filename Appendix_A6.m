%-------------------------------------------------------------------------
init;
show_current_script_name(mfilename('fullpath'));
%-------------------------------------------------------------------------

%% A6. Slope analysis between the inferred P(win) and the actual values: Regression analysis of slope across Experiments 1-4 and task modalities (LE, ES, and SP)

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

stat_A6_1 = []; % LE
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

cd1 = raw_Treatment == 1;
cd2 = raw_Treatment == 2;
cd3 = raw_Treatment == 3;
ce1 = raw_exp_num   == 1;
ce2 = raw_exp_num   == 2;
ce3 = raw_exp_num   == 3;
ce4 = raw_exp_num   == 4;

C1 = ce1;
C2 = (ce2);
C3 = (ce3);
C4 = (ce4);
C5 = (cd1);
C6 = (cd2);
C7 = (cd3);
C8 = ((ce2))&((cd2));
C9 = ((ce3))&((cd2));
C10 = ((ce4))&((cd2));
C11= ((ce2))&((cd3));
C12= ((ce3))&((cd3));
C13= ((ce4))&((cd3));

varNames = ["slope",'subject','T1','T2','T3','T4','M1','M2','M3'];
tbl = table(raw_slope,raw_subject,C1,C2,C3,C4,C5,C6,C7,'VariableNames',varNames);

formula = 'slope ~ 1 + (T2 + T3 + T4)*(M2 + M3) + (1 + (T2 + T3 + T4)*(M2 + M3)|subject)';
mdl = fitglme(tbl,formula);

tmp_CI = coefCI(mdl,'alpha',0.025);

tmp_coef = mdl.Coefficients.Estimate(:);
tmp_se   = mdl.Coefficients.SE(:);
tmp_tval = mdl.Coefficients.tStat(:);
tmp_pval = mdl.Coefficients.pValue(:);

tmp_CI98 = tmp_CI(:,:);
tmp_name = (mdl.Coefficients.Name);

for ii=1:length(tmp_name)
    stat_A6_1 = [stat_A6_1; tmp_name(ii) tmp_coef(ii) tmp_se(ii) tmp_tval(ii) tmp_pval(ii) tmp_CI98(ii,1) tmp_CI98(ii,2)];
end

%-------------------------------------------------------------------------

stat_A6_2 = []; %ES
formula = 'slope ~ 1 + (T2 + T3 + T4)*(M1 + M3) + (1 + (T2 + T3 + T4)*(M1 + M3)|subject)';
mdl = fitglme(tbl,formula);

tmp_CI = coefCI(mdl,'alpha',0.025);

tmp_coef = mdl.Coefficients.Estimate(:);
tmp_se   = mdl.Coefficients.SE(:);
tmp_tval = mdl.Coefficients.tStat(:);
tmp_pval = mdl.Coefficients.pValue(:);

tmp_CI98 = tmp_CI(:,:);
tmp_name = (mdl.Coefficients.Name);

for ii=1:length(tmp_name)
    stat_A6_2 = [stat_A6_2; tmp_name(ii) tmp_coef(ii) tmp_se(ii) tmp_tval(ii) tmp_pval(ii) tmp_CI98(ii,1) tmp_CI98(ii,2)];
end

%-------------------------------------------------------------------------

stat_A6_3 = []; %SP
formula = 'slope ~ 1 + (T2 + T3 + T4)*(M1 + M2) + (1 + (T2 + T3 + T4)*(M1 + M2)|subject)';
mdl = fitglme(tbl,formula);

tmp_CI = coefCI(mdl,'alpha',0.025);

tmp_coef = mdl.Coefficients.Estimate(:);
tmp_se   = mdl.Coefficients.SE(:);
tmp_tval = mdl.Coefficients.tStat(:);
tmp_pval = mdl.Coefficients.pValue(:);

tmp_CI98 = tmp_CI(:,:);
tmp_name = (mdl.Coefficients.Name);

for ii=1:length(tmp_name)
    stat_A6_3 = [stat_A6_3; tmp_name(ii) tmp_coef(ii) tmp_se(ii) tmp_tval(ii) tmp_pval(ii) tmp_CI98(ii,1) tmp_CI98(ii,2)];
end

%-------------------------------------------------------------------------

stat_A6_1
stat_A6_2
stat_A6_3