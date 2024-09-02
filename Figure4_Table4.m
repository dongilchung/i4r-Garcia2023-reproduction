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

%% Load stanDNCE Code 

stan_LE.alpha1 = [];
stan_LE.beta1  = [];
selected_exp = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
for ii=1:length(selected_exp)
    if ii>=6
        filenames = sprintf("data/fit/Stan/posterior_LE_models_5000permute_exp_%0.1f.mat",selected_exp(ii));
    else
        filenames = sprintf("data/fit/Stan/posterior_LE_models_5000permute_exp_%d.mat",selected_exp(ii));
    end
    tmp = load(filenames);
    stan_LE.alpha1 = [stan_LE.alpha1; tmp.medianPosteriorAlpha1];
    stan_LE.beta1  = [stan_LE.beta1; tmp.medianPosteriorBeta1];
end

stan_ES.beta1  = [];
stan_ES.midpoints1  = [];
selected_exp = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
for ii=1:length(selected_exp)
    if ii>=6
        filenames = sprintf("data/fit/Stan/new_posterior_ES_models_5000permute_exp_%0.1f.mat",selected_exp(ii));
    else
        filenames = sprintf("data/fit/Stan/new_posterior_ES_models_5000permute_exp_%d.mat",selected_exp(ii));
    end
    tmp = load(filenames);
    stan_ES.beta1  = [stan_ES.beta1; tmp.medianPosteriorBeta1];
    medianPosteriorMidPoints1 = squeeze(median(tmp.rawPosteriorMidPoints1,1));

    tmp2= zeros(size(medianPosteriorMidPoints1,1),8);
    tmp2(:,1:size(medianPosteriorMidPoints1,2)) = medianPosteriorMidPoints1;

    stan_ES.midpoints1  = [stan_ES.midpoints1; tmp2];
end

stan_EE.beta1  = [];
stan_EE.midpoints1  = [];
selected_exp = [5, 6.1, 6.2, 7.1 ,7.2, 8.1, 8.2, 9.1, 9.2];
for ii=1:length(selected_exp)
    if ii>=2
        filenames = sprintf("data/fit/Stan/new_posterior_EE_models_5000permute_exp_%0.1f.mat",selected_exp(ii));
    else
        filenames = sprintf("data/fit/Stan/new_posterior_EE_models_5000permute_exp_%d.mat",selected_exp(ii));
    end
    tmp = load(filenames);
    stan_EE.beta1  = [stan_EE.beta1; tmp.medianPosteriorBeta1];
    medianPosteriorMidPoints1 = squeeze(median(tmp.rawPosteriorMidPoints1,1));

    tmp2= zeros(size(medianPosteriorMidPoints1,1),8);
    tmp2(:,1:size(medianPosteriorMidPoints1,2)) = medianPosteriorMidPoints1;

    stan_EE.midpoints1  = [stan_EE.midpoints1; tmp2];
end

%% Figure4
%%% 1)LE
% Alpha 
figure,
plot(orig_data_LE_alpha, stan_LE.alpha1,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_LE_alpha, stan_LE.alpha1,'type','pearson')

p<10^(-16)

% Beta
figure,
plot(orig_data_LE_beta, stan_LE.beta1,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_LE_beta, stan_LE.beta1,'type','pearson')

p<10^(-16)

%%% 2)ES
% Beta
figure,
plot(orig_data_ES_beta, stan_ES.beta1,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_ES_beta, stan_ES.beta1,'type','pearson')

p<10^(-16)

% Midpoints
figure,
plot(orig_data_ES_midpoints(:), stan_ES.midpoints1(:),'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_ES_midpoints(:), stan_ES.midpoints1(:),'type','pearson')

p<10^(-16)

Table4_1 = [];
for ii=1:8
    figure,
    plot(orig_data_ES_midpoints(:,ii), stan_ES.midpoints1(:,ii),'k.','markersize',15)
    box off, set(gca,'LineWidth',2.0)
    [r,p] = corr(orig_data_ES_midpoints(:,ii), stan_ES.midpoints1(:,ii),'type','pearson')
    Table4_1 = [Table4_1; r p];
end

%%% 3)EE
% Beta
figure,
plot(orig_data_EE_beta, stan_EE.beta1,'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_EE_beta, stan_EE.beta1,'type','pearson')

p<10^(-16)

% Midpoints
figure,
plot(orig_data_EE_midpoints(:), stan_EE.midpoints1(:),'k.','markersize',15)
box off, set(gca,'LineWidth',2.0) 
[r,p] = corr(orig_data_EE_midpoints(:), stan_EE.midpoints1(:),'type','pearson')

p<10^(-2)

Table4_2 = [];
for ii=1:8
    figure,
    plot(orig_data_EE_midpoints(:,ii), stan_EE.midpoints1(:,ii),'k.','markersize',15)
    box off, set(gca,'LineWidth',2.0)
    [r,p] = corr(orig_data_EE_midpoints(:,ii), stan_EE.midpoints1(:,ii),'type','pearson')
    Table4_2 = [Table4_2; r p];
end
