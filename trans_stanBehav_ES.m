clc
clear all
close all

%%
init;
selected_exp = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];

stanAccum_choice   = [];
stanAccum_subject  = [];
stanAccum_ntrial   = [];
stanAccum_p1_sym   = [];
stanAccum_p2_lot   = [];
stanAccum_alt_prob = [];
stanAccum_nsub     = [];
stanAccum_ncon     = [];
stanAccum_session  = [];
stanAccum_exp_num  = [];
for exp_num = selected_exp
    data = de.extract_ES(exp_num);
    fit_params.cho = data.cho;
    fit_params.p1  = data.p1; % symbol
    fit_params.p2  = data.p2; % lottery
    fit_params.ntrials = size(data.cho, 2);
    fit_params.nsub = data.nsub;
    fit_params.sess = data.sess;
    fit_params.exp_num = num2str(exp_num);

    fit_params.subject = [];
    for ii=1:size(fit_params.cho,1)
        fit_params.subject = [fit_params.subject; ii.*ones(1,size(fit_params.cho,2))];
    end

    stanAccum_choice   = [stanAccum_choice; reshape(fit_params.cho',[size(fit_params.cho,1).*size(fit_params.cho,2),1])];
    stanAccum_p1_sym   = [stanAccum_p1_sym; reshape(fit_params.p1',[size(fit_params.cho,1).*size(fit_params.cho,2),1])];
    stanAccum_p2_lot   = [stanAccum_p2_lot; reshape(fit_params.p2',[size(fit_params.cho,1).*size(fit_params.cho,2),1])];
    stanAccum_subject  = [stanAccum_subject; reshape(fit_params.subject',[size(fit_params.subject,1).*size(fit_params.subject,2),1])];    
    stanAccum_ntrial   = [stanAccum_ntrial; fit_params.ntrials.*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    stanAccum_nsub     = [stanAccum_nsub; fit_params.nsub.*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    stanAccum_session  = [stanAccum_session; fit_params.sess.*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    if isnumeric(fit_params.exp_num)
        stanAccum_exp_num  = [stanAccum_exp_num; (fit_params.exp_num).*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    else
        stanAccum_exp_num  = [stanAccum_exp_num; str2num(fit_params.exp_num).*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    end
end

save('stan/Data/stanBehav_ES_choices.mat','stanAccum_*')