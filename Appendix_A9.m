%-------------------------------------------------------------------------
init;
show_current_script_name(mfilename('fullpath'));
%-------------------------------------------------------------------------

%% A9. Slope analysis between the inferred P(win) and the actual values: T-test comparison of slopes across Experiments 6.1-6.2 and task modalities (ES, EE)

%-------------------------------------------------------------------------%
% parameters of the script                                                %
%-------------------------------------------------------------------------%
selected_exp = [5, 6.1, 6.2];
modalities = {'ES', 'EE'};
displayfig = 'on';
colors = [orange;green];
% filenames
filename = 'Fig3D';
figfolder = 'fig';

figname = sprintf('%s/%s.svg', figfolder, filename);
stats_filename = sprintf('data/stats/%s.csv', filename);

fitname = 'data/fit/midpoints_%s_%s_session_%d';

%-------------------------------------------------------------------------%
% prepare data                                                            %
%-------------------------------------------------------------------------%
% stats_data is table that is used to compute stats later
stats_data = table();

num = 0;
sub_count = 0;
for exp_num = selected_exp
    num = num + 1;
    
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
                [midpoints(mod_num, :, :), throw] = get_qvalues(sim_params);
                
            case {'EE', 'ES'}
                param = load(...
                    sprintf(fitname,...
                    modalities{mod_num}, name , sess));
                
                midpoints(mod_num, :, :) = param.midpoints;
                
            case 'PM'
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
                sub+sub_count, exp_num, slope(mod_num, sub, 2),...
                {modalities{mod_num}}, 'variablenames',...
                {'subject', 'exp_num', 'slope', 'modality'}...
                );
            stats_data = [stats_data; T1];
        end
    end
    sub_count = sub_count+sub;    
end

%-------------------------------------------------------------------------

raw_slope   = stats_data.slope(stats_data.exp_num~=5);
raw_exp_num = stats_data.exp_num(stats_data.exp_num~=5);
raw_subject = stats_data.subject(stats_data.exp_num~=5);
tmp = stats_data.modality(stats_data.exp_num~=5);
raw_Treatment = [];
for ii=1:length(raw_slope)
    switch tmp{ii}
        case 'ES'
            raw_Treatment = [raw_Treatment; 0];
        case 'EE'
            raw_Treatment = [raw_Treatment; 1];
    end
end
stat_A9_1 = []; %ES

cd1 = raw_Treatment == 0;
cd2 = raw_Treatment == 1;
ce1 = raw_exp_num   == 6.1;
ce2 = raw_exp_num   == 6.2;

C1 = ce1;
C2 = (ce2);
C3 = (cd1);
C4 = (cd2);

varNames = ["slope",'subject','T1','T2','M1','M2'];
tbl = table(raw_slope,raw_subject,C1,C2,C3,C4,'VariableNames',varNames);

formula = 'slope ~ 1 + (T2)*(M2) + (1 + (T2)*(M2)|subject)';
mdl = fitglme(tbl,formula);

tmp_CI = coefCI(mdl,'alpha',0.025);

tmp_coef = mdl.Coefficients.Estimate(:);
tmp_se   = mdl.Coefficients.SE(:);
tmp_tval = mdl.Coefficients.tStat(:);
tmp_pval = mdl.Coefficients.pValue(:);

tmp_CI98 = tmp_CI(:,:);
tmp_name = (mdl.Coefficients.Name);

for ii=1:length(tmp_name)
    stat_A9_1 = [stat_A9_1; tmp_name(ii) tmp_coef(ii) tmp_se(ii) tmp_tval(ii) tmp_pval(ii) tmp_CI98(ii,1) tmp_CI98(ii,2)];
end

stat_A9_2 = []; %EE

cd1 = raw_Treatment == 0;
cd2 = raw_Treatment == 1;
ce1 = raw_exp_num   == 6.1;
ce2 = raw_exp_num   == 6.2;

C1 = ce1;
C2 = (ce2);
C3 = (cd1);
C4 = (cd2);

varNames = ["slope",'subject','T1','T2','M1','M2'];
tbl = table(raw_slope,raw_subject,C1,C2,C3,C4,'VariableNames',varNames);

formula = 'slope ~ 1 + (T2)*(M1) + (1 + (T2)*(M1)|subject)';
mdl = fitglme(tbl,formula);

tmp_CI = coefCI(mdl,'alpha',0.025);

tmp_coef = mdl.Coefficients.Estimate(:);
tmp_se   = mdl.Coefficients.SE(:);
tmp_tval = mdl.Coefficients.tStat(:);
tmp_pval = mdl.Coefficients.pValue(:);

tmp_CI98 = tmp_CI(:,:);
tmp_name = (mdl.Coefficients.Name);

for ii=1:length(tmp_name)
    stat_A9_2 = [stat_A9_2; tmp_name(ii) tmp_coef(ii) tmp_se(ii) tmp_tval(ii) tmp_pval(ii) tmp_CI98(ii,1) tmp_CI98(ii,2)];
end

%-------------------------------------------------------------------------

stat_A9_1
stat_A9_2