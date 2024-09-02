clc
clear all
close all

%%
init;
selected_exp = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
sessions = [0, 1];

stanAccum_choice   = [];
stanAccum_subject  = [];
stanAccum_ntrial   = [];
stanAccum_con      = [];
stanAccum_cf       = [];
stanAccum_out      = [];
stanAccum_cfout    = [];
stanAccum_nsub     = [];
stanAccum_ncon     = [];
stanAccum_session  = [];
stanAccum_exp_num  = [];
for exp_num = selected_exp
    data = de.extract_LE(exp_num);
    fit_params.cho = data.cho;
    fit_params.cfcho = data.cfcho; % just reversed cho
    fit_params.out = data.out==1;
    fit_params.cfout = data.cfout==1;
    fit_params.con = data.con;
    fit_params.fit_cf = (exp_num>2);
    fit_params.ntrials = size(data.cho, 2);
    fit_params.model = 1;
    fit_params.nsub = data.nsub;
    fit_params.sess = data.sess;
    fit_params.exp_num = num2str(exp_num);
    fit_params.decision_rule = 1;
    fit_params.q = 0.5;
    fit_params.noptions = 2;
    fit_params.ncond = length(unique(data.con));

    fit_params.subject = [];
    for ii=1:size(fit_params.cho,1)
        fit_params.subject = [fit_params.subject; ii.*ones(1,size(fit_params.cho,2))];
    end
    stanAccum_choice   = [stanAccum_choice; reshape(fit_params.cho',[size(fit_params.cho,1).*size(fit_params.cho,2),1])];
    stanAccum_subject  = [stanAccum_subject; reshape(fit_params.subject',[size(fit_params.subject,1).*size(fit_params.subject,2),1])];    
    stanAccum_ntrial   = [stanAccum_ntrial; fit_params.ntrials.*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    stanAccum_con      = [stanAccum_con; reshape(fit_params.con',[size(fit_params.cho,1).*size(fit_params.cho,2),1])];
    stanAccum_cf       = [stanAccum_cf; fit_params.fit_cf.*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    stanAccum_out      = [stanAccum_out; reshape(fit_params.out',[size(fit_params.cho,1).*size(fit_params.cho,2),1])];
    stanAccum_cfout    = [stanAccum_cfout; reshape(fit_params.cfout',[size(fit_params.cho,1).*size(fit_params.cho,2),1])];
    stanAccum_nsub     = [stanAccum_nsub; fit_params.nsub.*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    stanAccum_ncon     = [stanAccum_ncon; fit_params.ncond.*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    stanAccum_session  = [stanAccum_session; fit_params.sess.*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    if isnumeric(fit_params.exp_num)
        stanAccum_exp_num  = [stanAccum_exp_num; (fit_params.exp_num).*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    else
        stanAccum_exp_num  = [stanAccum_exp_num; str2num(fit_params.exp_num).*ones(size(fit_params.cho,1).*size(fit_params.cho,2),1)];
    end
end

save('stan/Data/stanBehav_LE_choices.mat','stanAccum_*')