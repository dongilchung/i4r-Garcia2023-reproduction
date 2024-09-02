%-------------------------------------------------------------------------%
clc
clear all
close all

init;
show_current_script_name(mfilename('fullpath'));

%-------------------------------------------------------------------------%
%% Load Garcia Data
selected_exp_LE = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
fitname_LE = 'data/fit/Original_Paper/learning_LE_%s_session_%d';
orig_data_LE_alpha = [];
orig_data_LE_beta  = [];
for exp_num = selected_exp_LE
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    nsub = de.get_nsub_from_exp_num(exp_num);
    
    throw = de.extract_ES(exp_num);
    nsym = length(unique(throw.p1));
    p1 = unique(throw.p1)'.*100;
 
    
    sim_params.exp_num = exp_num;
    sim_params.de = de;
    sim_params.sess = sess;
    sim_params.exp_name = name;
    sim_params.nsub = nsub;
    sim_params.path = fitname_LE;
                    
    sim_params.model = 1;
    [midpoints1, throw] = get_qvalues(sim_params);
    orig_data_LE_alpha = [orig_data_LE_alpha; throw.alpha1];
    orig_data_LE_beta  = [orig_data_LE_beta; throw.beta1];
end
orig_data_ES_beta = [];
selected_exp_ES = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
fitname_ES = 'data/fit/Original_Paper/midpoints_ES_%s_session_%d';
for exp_num = selected_exp_ES
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_ES, name, sess));
    orig_data_ES_beta  = [orig_data_ES_beta; param.beta1];
end
orig_data_ES_midpoints = [];
selected_exp_ES = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
fitname_ES = 'data/fit/Original_Paper/midpoints_ES_%s_session_%d';
for exp_num = selected_exp_ES
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_ES, name, sess));
    tmp = zeros(size(param.midpoints,1),8);
    tmp(:,1:size(param.midpoints,2)) = param.midpoints;
    orig_data_ES_midpoints  = [orig_data_ES_midpoints; tmp];
end
orig_data_EE_beta = [];
selected_exp_EE = [5, 6.1, 6.2, 7.1 ,7.2, 8.1, 8.2, 9.1, 9.2];
fitname_EE = 'data/fit/Original_Paper/midpoints_EE_%s_session_%d';
for exp_num = selected_exp_EE
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_EE, name, sess));
    orig_data_EE_beta  = [orig_data_EE_beta; param.beta1];
end
orig_data_EE_midpoints = [];
selected_exp_EE = [5, 6.1, 6.2, 7.1 ,7.2, 8.1, 8.2, 9.1, 9.2];
fitname_EE = 'data/fit/Original_Paper/midpoints_EE_%s_session_%d';
for exp_num = selected_exp_EE
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_EE, name, sess));
    tmp = zeros(size(param.midpoints,1),8);
    tmp(:,1:size(param.midpoints,2)) = param.midpoints;
    orig_data_EE_midpoints  = [orig_data_EE_midpoints; tmp];
end

%% Load DNCE Code

repl_LE = load('data/repl_parameter_recovery/mh_fit_LE');
repl_EE = load('data/repl_parameter_recovery/mh_fit_EE');
repl_ES = load('data/repl_parameter_recovery/mh_fit_ES');

%% Figure 3
%%% 1)LE
% Alpha 
figure,
plot(orig_data_LE_alpha, repl_LE.repl_alpha,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_LE_alpha, repl_LE.repl_alpha,'type','pearson')

p<10^(-16)

% Beta
figure,
plot(orig_data_LE_beta, repl_LE.repl_beta,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_LE_beta, repl_LE.repl_beta,'type','pearson')

p<10^(-16)

%%% 2)ES
% Beta
figure,
plot(orig_data_ES_beta, repl_ES.repl_beta,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_ES_beta, repl_ES.repl_beta,'type','pearson')

p<10^(-16)

% Midpoints
figure,
plot(orig_data_ES_midpoints(:), repl_ES.repl_midPoint(:),'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_ES_midpoints(:), repl_ES.repl_midPoint(:),'type','pearson')

p<10^(-16)

Table3_1 =[];
for ii=1:8
    figure,
    plot(orig_data_ES_midpoints(:,ii), repl_ES.repl_midPoint(:,ii),'k.','markersize',15)
    box off, set(gca,'LineWidth',2.0)
    [r,p] = corr(orig_data_ES_midpoints(:,ii), repl_ES.repl_midPoint(:,ii),'type','pearson')
    Table3_1 = [Table3_1; r p];
end

%%% 3)EE
% Beta
figure,
plot(orig_data_EE_beta, repl_EE.repl_beta,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_EE_beta, repl_EE.repl_beta,'type','pearson')

% Midpoints
figure,
plot(orig_data_EE_midpoints(:), repl_EE.repl_midPoint(:),'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_EE_midpoints(:), repl_EE.repl_midPoint(:),'type','pearson')

Table3_2 =[];
for ii=1:8
    figure,
    plot(orig_data_EE_midpoints(:,ii), repl_EE.repl_midPoint(:,ii),'k.','markersize',15)
    box off, set(gca,'LineWidth',2.0)
    [r,p] = corr(orig_data_EE_midpoints(:,ii), repl_EE.repl_midPoint(:,ii),'type','pearson')
    Table3_2 = [Table3_2; r p];
end

