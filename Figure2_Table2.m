%-------------------------------------------------------------------------%
clc
clear all
close all

init;
show_current_script_name(mfilename('fullpath'));

%-------------------------------------------------------------------------%
%% Load Garcia Code 
selected_exp_LE = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
fitname_LE = 'data/fit/learning_LE_%s_session_%d';
orig_code_LE_alpha = [];
orig_code_LE_beta  = [];
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
    orig_code_LE_alpha = [orig_code_LE_alpha; throw.alpha1];
    orig_code_LE_beta  = [orig_code_LE_beta; throw.beta1];
end
orig_code_ES_beta = [];
selected_exp_ES = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
fitname_ES = 'data/fit/midpoints_ES_%s_session_%d';
for exp_num = selected_exp_ES
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_ES, name, sess));
    orig_code_ES_beta  = [orig_code_ES_beta; param.beta1];
end
orig_code_ES_midpoints = [];
selected_exp_ES = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
fitname_ES = 'data/fit/midpoints_ES_%s_session_%d';
for exp_num = selected_exp_ES
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_ES, name, sess));
    tmp = zeros(size(param.midpoints,1),8);
    tmp(:,1:size(param.midpoints,2)) = param.midpoints;
    orig_code_ES_midpoints  = [orig_code_ES_midpoints; tmp];
end
orig_code_EE_beta = [];
selected_exp_EE = [5, 6.1, 6.2, 7.1 ,7.2, 8.1, 8.2, 9.1, 9.2];
fitname_EE = 'data/fit/midpoints_EE_%s_session_%d';
for exp_num = selected_exp_EE
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_EE, name, sess));
    orig_code_EE_beta  = [orig_code_EE_beta; param.beta1];
end
orig_code_EE_midpoints = [];
selected_exp_EE = [5, 6.1, 6.2, 7.1 ,7.2, 8.1, 8.2, 9.1, 9.2];
fitname_EE = 'data/fit/midpoints_EE_%s_session_%d';
for exp_num = selected_exp_EE
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_EE, name, sess));
    tmp = zeros(size(param.midpoints,1),8);
    tmp(:,1:size(param.midpoints,2)) = param.midpoints;
    orig_code_EE_midpoints  = [orig_code_EE_midpoints; tmp];
end

%% Load Garcia Code - Parameter recovery
selected_exp_LE = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
fitname_LE = 'data/pr/learning_LE%s_session_%d';
orig_pr_code_LE_alpha = [];
orig_pr_code_LE_beta  = [];
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
    orig_pr_code_LE_alpha = [orig_pr_code_LE_alpha; throw.alpha1];
    orig_pr_code_LE_beta  = [orig_pr_code_LE_beta; throw.beta1];
end
orig_pr_code_ES_beta = [];
selected_exp_ES = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
fitname_ES = 'data/pr/midpoints_ES_%s_session_%d';
for exp_num = selected_exp_ES
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_ES, name, sess));
    orig_pr_code_ES_beta  = [orig_pr_code_ES_beta; param.beta1];
end
orig_pr_code_ES_midpoints = [];
selected_exp_ES = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
fitname_ES = 'data/pr/midpoints_ES_%s_session_%d';
for exp_num = selected_exp_ES
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_ES, name, sess));
    tmp = zeros(size(param.midpoints,1),8);
    tmp(:,1:size(param.midpoints,2)) = param.midpoints;
    orig_pr_code_ES_midpoints  = [orig_pr_code_ES_midpoints; tmp];
end
orig_pr_code_EE_beta = [];
selected_exp_EE = [5, 6.1, 6.2, 7.1 ,7.2, 8.1, 8.2, 9.1, 9.2];
fitname_EE = 'data/pr/midpoints_EE_%s_session_%d';
for exp_num = selected_exp_EE
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_EE, name, sess));
    orig_pr_code_EE_beta  = [orig_pr_code_EE_beta; param.beta1];
end
orig_pr_code_EE_midpoints = [];
selected_exp_EE = [5, 6.1, 6.2, 7.1 ,7.2, 8.1, 8.2, 9.1, 9.2];
fitname_EE = 'data/pr/midpoints_EE_%s_session_%d';
for exp_num = selected_exp_EE
    sess = de.get_sess_from_exp_num(exp_num);
    name = de.get_name_from_exp_num(exp_num);
    param = load(sprintf(fitname_EE, name, sess));
    tmp = zeros(size(param.midpoints,1),8);
    tmp(:,1:size(param.midpoints,2)) = param.midpoints;
    orig_pr_code_EE_midpoints  = [orig_pr_code_EE_midpoints; tmp];
end

%% Figure 2
%%% 1)LE
% Alpha 
figure,
plot(orig_code_LE_alpha, orig_pr_code_LE_alpha,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_code_LE_alpha, orig_pr_code_LE_alpha,'type','pearson')

p<10^(-16)

% Beta
figure,
plot(orig_code_LE_beta, orig_pr_code_LE_beta,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_code_LE_beta, orig_pr_code_LE_beta,'type','pearson')

p<10^(-16)

%%% 2)ES
% Beta
figure,
plot(orig_code_ES_beta, orig_pr_code_ES_beta,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_code_ES_beta, orig_pr_code_ES_beta,'type','pearson')

p<10^(-16)

% Midpoints
figure,
plot(orig_code_ES_midpoints(:), orig_pr_code_ES_midpoints(:),'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_code_ES_midpoints(:), orig_pr_code_ES_midpoints(:),'type','pearson')

p<10^(-16)

Table2_1 = [];
for ii=1:8
    figure,
    plot(orig_code_ES_midpoints(:,ii), orig_pr_code_ES_midpoints(:,ii),'k.','markersize',15)
    box off, set(gca,'LineWidth',2.0)
    [r,p] = corr(orig_code_ES_midpoints(:,ii), orig_pr_code_ES_midpoints(:,ii),'type','pearson')
    Table2_1 = [Table2_1; r p];
end

%%% 3)EE
% Beta
figure,
plot(orig_code_EE_beta, orig_pr_code_EE_beta,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_code_EE_beta, orig_pr_code_EE_beta,'type','pearson')

p<10^(-16)

% Midpoints
figure,
plot(orig_code_EE_midpoints(:), orig_pr_code_EE_midpoints(:),'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_code_EE_midpoints(:), orig_pr_code_EE_midpoints(:),'type','pearson')

p<10^(-16)

Table2_2 = [];
for ii=1:8
    figure,
    plot(orig_code_EE_midpoints(:,ii), orig_pr_code_EE_midpoints(:,ii),'k.','markersize',15)
    box off, set(gca,'LineWidth',2.0)
    [r,p] = corr(orig_code_EE_midpoints(:,ii), orig_pr_code_EE_midpoints(:,ii),'type','pearson')
    Table2_2 = [Table2_2; r p];
end
