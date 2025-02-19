classdef DataExtractionCSV < handle

    properties (SetAccess = public)
        d
        filenames
    end

    methods (Static)

        function [data, sub_ids] = get_data(filename)

            data = readtable(filename);
            sub_ids = unique(data.sub_id);
            sub_ids = sub_ids(~isnan(sub_ids));

        end
            function error_exclude = extract_learning_data(data, sub_ids, session)
            % check that needed fields have same length at init time
            i = 1;
            error_exclude = [];

            for id = 1:length(sub_ids)
                try

                    
                    sub = sub_ids(id);
                    mask_sub = data.sub_id  == sub;
                    mask_sess = ismember(data.sess, session);
                    %mask_eli = data(:).elic) == -1;
                    mask = logical(mask_sub .* mask_sess);

                    [noneed, trialorder] = sort(data(mask,:).trial);

                    tempcho = data(mask,:).cho;
                    cho(i, :) = tempcho(trialorder);

                    cfcho(i, :) = 3 - cho(i, :);

                    tempout = data(mask,:).out;
                    out(i, :) = tempout(trialorder);
                    tempcorr = data(mask,:).corr;

                    corr(i, :) = tempcorr(trialorder);
                    temprew = data(mask,:).rew;

                    rew(i, :) = temprew(trialorder);

                    tempcon = data(mask,:).cond;
                    con(i, :) = tempcon(trialorder) + 1;

                    tempcfout = data(mask,:).cfout;
                    cfout(i, :) = tempcfout(trialorder);

                    temp_p1 = data(mask,:).p1;
                    p1(i, :) = temp_p1(trialorder);

                    temp_p2 = data(mask,:).p2;
                    p2(i, :) = temp_p2(trialorder);

                    temp_rtime = data(mask,:).rtime;
                    rtime(i, :) = temp_rtime(trialorder);

                    temp_ev1 = data(mask,:).ev1;
                    ev1(i, :) = temp_ev1(trialorder);

                    temp_ev2 = data(mask,:).ev2;
                    ev2(i, :) = temp_ev2(trialorder);

                    i = i + 1;
                catch e
                    error_exclude(length(error_exclude) + 1) = i;

                    fprintf(1, '\n There has been an error while treating subject %d \n', i);
                    fprintf(1,'\n The identifier was: %s \n',e.identifier);
                    fprintf(1,'\n The message was: %s \n',e.message);
                end
            end
        end



        function to_keep = exclude_subjects(data, sub_ids,...
                ES_catch_trial_threshold, PM_catch_trial_threshold, PM_corr_threshold, ...
                rtime_threshold, allowed_nb_of_rows)

            to_keep = [];
            sums = [];
            i = 1;
            n_complete = 0;
            possible_eli = {'ES', 'EE', 'SP'};
            for id = 1:length(sub_ids)
                sub = sub_ids(id);
                sums(id) = sum(data.sub_id == sub);

                if ismember(sum(data.sub_id == sub), allowed_nb_of_rows)
                    n_complete = n_complete + 1;
                    for eli = 1:length(possible_eli)

                        % if ES
                        if (strcmp(possible_eli{eli}, 'ES')) || (strcmp(possible_eli{eli}, 'EE'))
                            mask_eli = strcmp(data.phase, possible_eli{eli});

                            mask_sub = data.sub_id == sub;
                            mask_catch_trial = data.catch_trial == 1;
                            mask_sess = ismember(data.sess, [0,1]);
                            mask = logical(mask_sub .* mask_sess .* mask_catch_trial .* mask_eli);

                            temp_corr = data(mask,:).corr;

                            corr_catch_trial{i, eli} = temp_corr;

                            mask = logical(mask_sub .* mask_sess .* mask_eli);
                            rtime{i, eli} = data(mask,:).rtime;
                        
                        % if SP
                        else
                            mask_sub = data.sub_id == sub;
                            mask_sess = ismember(data.sess, [0,1]);
                            mask_eli = strcmp(data.phase, 'SP');
                            mask_catch = (data.catch_trial == 1) | (data.catch_trial == 0);
                            mask = logical(mask_sub .* mask_sess .* mask_eli .* mask_catch);
                            rtime{i, eli} = data(mask,:).rtime;
                        end
                    end
                    
                    if (mean(corr_catch_trial{i, 1}) >= ES_catch_trial_threshold) &&...%(mean(corr_catch_trial{i, 2}) >= PM_catch_trial_threshold)...                     
                            (sum(vertcat(rtime{i, :}) > rtime_threshold) < 1) % && (sum(corr1{i, 3}) > 0)
                        to_keep(length(to_keep) + 1) = sub;
                    
                    end
                    i = i + 1;
                    
                end

            end
            fprintf('N = %d \n', length(sub_ids)); 
            fprintf('N complete = %d \n', n_complete); 
            fprintf('N after exclusion = %d \n', length(to_keep)); 


        end

    end


    methods

        function obj = DataExtractionCSV(d, filenames)
            obj.d = d;
            obj.filenames = filenames;
        end

        function name = get_name_from_exp_num(obj, exp_num)
            % load data
            name = char(obj.filenames{round(exp_num(1))});
        end

        function nsub = get_nsub_from_exp_num(obj, exp_num)
            nsub = obj.d.(get_name_from_exp_num(obj, exp_num)).nsub;
        end

        function sess = get_sess_from_exp_num(obj, exp_num)
            if length(exp_num) == 1
                if (length(num2str(exp_num)) == 1) && (exp_num >= 6)
                    sess = [0,1];
                else

                    sess = round((exp_num - round(exp_num)) * 10 - 1);
                    sess = sess .* (sess ~= -1);
                end

            elseif length(exp_num) == 2
                sess = [0,1];

            else
                error('Problem matching session');
            end

        end


        function finalstruct = merge2sess(obj, structlist)
            f = fieldnames(structlist(1));
            for i = 1:length(f)
                if ~ismember(f{i}, {'sess', 'nsub', 'name', 'exp_num'})
                    structlist(1).(f{i}) = [structlist(1).(f{i}) structlist(2).(f{i})];
                end
            end
            finalstruct = structlist(1);
        end

        function [data, sub_ids, sess, name, nsub] = prepare(obj, exp_num)
            sess = obj.get_sess_from_exp_num(exp_num);
            name = obj.get_name_from_exp_num(exp_num);
            nsub = obj.get_nsub_from_exp_num(exp_num);
            data = obj.d.(name).data;
            sub_ids = obj.d.(name).sub_ids;
        end

        function new_data = extract_LE(obj, exp_num)
            [data, sub_ids, session, name,nsub] = prepare(obj, exp_num);
            new_data = struct();

            new_data.sess = session;
            new_data.name = name;
            new_data.nsub = nsub;
            new_data.exp_num = exp_num;
            new_data.sub_ids = sub_ids;


            for isess = 1:length(session)
                i = 1;

                for id = 1:length(sub_ids)
                    try
                        sub = sub_ids(id);
                        
                        mask_eli = strcmp(data.phase, 'LE');
                        mask_sub = data.sub_id == sub;
                        mask_catch_trial = data.catch_trial == 0;

                        mask_sess = data.sess ==  session(isess);
                        mask = logical(mask_sub .* mask_sess .* mask_eli .* mask_catch_trial);

                        trialorder = data(mask,:).trial;

                        if ~issorted(trialorder)
                            [noneed, trialorder] = sort(data(mask,:).trial);
                        else
                            trialorder = 1:length(trialorder);
                        end


                        temp_corr = data(mask, :).corr;
                        new_data.corr(i, :) = temp_corr(trialorder);

                         temp_corr = data(mask, :).cond+1;
                        new_data.con(i, :) = temp_corr(trialorder);

                        temp_cho = data(mask, :).cho;
                        new_data.cho(i, :) = temp_cho(trialorder);

                        new_data.cfcho(i, :) = 3 - new_data.cho(i, :);

                        temp_out = data(mask, :).out;
                        new_data.out(i, :) = temp_out(trialorder);

                        temp_cfout = data(mask, :).cfout;
                        new_data.cfout(i, :) = temp_cfout(trialorder);

                        temp_ev1 = data(mask, :).ev1;                       

                        new_data.ev1(i, :) = temp_ev1(trialorder);

                        temp_catch_trial = data(mask, :).catch_trial;
                        new_data.ctch(i, :) = temp_catch_trial(trialorder);

%                         new_data.ctch_p1(i, :) = data(mask2, obj.idx.p1);
% 
%                         new_data.ctch_p2(i, :) = data(mask2, obj.idx.p2);
% 
%                         new_data.ctch_corr(i, :) = data(mask2, obj.idx.corr);

%                         temp_cont1 = data(mask, obj.idx.cont1);
%                         new_data.cont1(i, :) = temp_cont1(trialorder);

                        temp_ev2 = data(mask,:).ev2;
                        new_data.ev2(i, :) = temp_ev2(trialorder);

%                         temp_cont2 = data(mask, obj.idx.cont2);
%                         new_data.cont2(i, :) = temp_cont2(trialorder);

                        temp_p1 = data(mask,:).p1;
                        new_data.p1(i, :) = temp_p1(trialorder);

                        temp_p2 = data(mask, :).p2;
                        new_data.p2(i, :) = temp_p2(trialorder);

%                         temp_dist = data(mask,:);
%                         new_data.dist(i, :) = temp_dist(trialorder)./100;

                        temp_rtime = data(mask, :).rtime;
                        new_data.rtime(i, :) = temp_rtime(trialorder);
% 
                        temp_trial = data(mask,:).trial;
                        new_data.trial(i, :) = temp_trial(trialorder);

%                         new_data.catch_trial(i, :) = data(mask2, obj.idx.catch_trial);                        


                        i = i + 1;

                    catch e
                        fprintf(1, '\n There has been an error while treating subject %d \n', i);
                        fprintf(1,'The identifier was:\n%s',e.identifier);
                        fprintf(1,'There was an error! The message was:\n%s',e.message);
                    end
                end

                structlist(isess) = new_data;

            end

            if length(session) > 1
                new_data = obj.merge2sess(structlist);
            end

        end


        function new_data = extract_ES(obj, exp_num)
            [data, sub_ids, session, name,nsub] = prepare(obj, exp_num);
            new_data = struct();

            new_data.sess = session;
            new_data.name = name;
            new_data.nsub = nsub;
            new_data.exp_num = exp_num;
            new_data.id = sub_ids;


            for isess = 1:length(session)
                i = 1;

                for id = 1:length(sub_ids)
                    try
                        sub = sub_ids(id);
                        
                        mask_eli = strcmp(data.phase, 'ES');
                        mask_sub = data.sub_id == sub;
                        mask_catch_trial = data.catch_trial == 0;
%                         mask_ycatch_trial = data(:, obj.idx.catch_trial) == 1;

                        % before exp. 5 op2 has value -1 in ED while after it
                        % takes value 0 (because there were no EE befor                   

                        mask_sess = data.sess ==  session(isess);
                        mask = logical(mask_sub .* mask_sess .* mask_eli .* mask_catch_trial);

                        
%                         mask2 = logical(mask_sub .* mask_sess .* mask_eli .* mask_vs_lot .* mask_ycatch_trial.* mask_vs_lot);

                        trialorder = data(mask,:).trial;

                        if ~issorted(trialorder)
                            [noneed, trialorder] = sort(data(mask,:).trial);
                        else
                            trialorder = 1:length(trialorder);
                        end


                        temp_corr = data(mask, :).corr;
                        new_data.corr(i, :) = temp_corr(trialorder);

                        temp_cho = data(mask, :).cho;
                        new_data.cho(i, :) = temp_cho(trialorder);

                        new_data.cfcho(i, :) = 3 - new_data.cho(i, :);

                        temp_out = data(mask, :).out;
                        new_data.out(i, :) = temp_out(trialorder);

                        temp_ev1 = data(mask, :).ev1;                       

                        new_data.ev1(i, :) = temp_ev1(trialorder);

                        temp_catch_trial = data(mask, :).catch_trial;
                        new_data.ctch(i, :) = temp_catch_trial(trialorder);

%                         new_data.ctch_p1(i, :) = data(mask2, obj.idx.p1);
% 
%                         new_data.ctch_p2(i, :) = data(mask2, obj.idx.p2);
% 
%                         new_data.ctch_corr(i, :) = data(mask2, obj.idx.corr);

%                         temp_cont1 = data(mask, obj.idx.cont1);
%                         new_data.cont1(i, :) = temp_cont1(trialorder);

                        temp_ev2 = data(mask,:).ev2;
                        new_data.ev2(i, :) = temp_ev2(trialorder);

%                         temp_cont2 = data(mask, obj.idx.cont2);
%                         new_data.cont2(i, :) = temp_cont2(trialorder);

                        temp_p1 = data(mask,:).p1;
                        new_data.p1(i, :) = temp_p1(trialorder);

                        temp_p2 = data(mask, :).p2;
                        new_data.p2(i, :) = temp_p2(trialorder);

%                         temp_dist = data(mask,:);
%                         new_data.dist(i, :) = temp_dist(trialorder)./100;

                        temp_rtime = data(mask, :).rtime;
                        new_data.rtime(i, :) = temp_rtime(trialorder);
% 
                        temp_trial = data(mask,:).trial;
                        new_data.trial(i, :) = temp_trial(trialorder);

%                         new_data.catch_trial(i, :) = data(mask2, obj.idx.catch_trial);                        

                        i = i + 1;

                    catch e
                        fprintf(1, '\n There has been an error while treating subject %d \n', i);
                        fprintf(1,'The identifier was:\n%s',e.identifier);
                        fprintf(1,'There was an error! The message was:\n%s',e.message);
                    end
                end

                structlist(isess) = new_data;

            end

            if length(session) > 1
                new_data = obj.merge2sess(structlist);
            end

        end

        function new_data = extract_EE(obj, exp_num)
            [data, sub_ids, session, name,nsub] = prepare(obj, exp_num);
            new_data = struct();

            new_data.sess = session;
            new_data.name = name;
            new_data.nsub = nsub;
            new_data.exp_num = exp_num;
            new_data.id = sub_ids;


            for isess = 1:length(session)
                i = 1;

                for id = 1:length(sub_ids)
                    try
                        sub = sub_ids(id);
                        
                        mask_eli = strcmp(data.phase, 'EE');
                        mask_sub = data.sub_id == sub;
                        mask_catch_trial = data.catch_trial == 0;
%                         mask_ycatch_trial = data(:, obj.idx.catch_trial) == 1;

                        % before exp. 5 op2 has value -1 in ED while after it
                        % takes value 0 (because there were no EE befor                   

                        mask_sess = data.sess ==  session(isess);
                        mask = logical(mask_sub .* mask_sess .* mask_eli .* mask_catch_trial);

                        
%                         mask2 = logical(mask_sub .* mask_sess .* mask_eli .* mask_vs_lot .* mask_ycatch_trial.* mask_vs_lot);

                        trialorder = data(mask,:).trial;

                        if ~issorted(trialorder)
                            [noneed, trialorder] = sort(data(mask,:).trial);
                        else
                            trialorder = 1:length(trialorder);
                        end


                        temp_corr = data(mask, :).corr;
                        new_data.corr(i, :) = temp_corr(trialorder);

                        temp_cho = data(mask, :).cho;
                        new_data.cho(i, :) = temp_cho(trialorder);

                        new_data.cfcho(i, :) = 3 - new_data.cho(i, :);

                        temp_out = data(mask, :).out;
                        new_data.out(i, :) = temp_out(trialorder);

                        temp_ev1 = data(mask, :).ev1;                       

                        new_data.ev1(i, :) = temp_ev1(trialorder);

                        temp_catch_trial = data(mask, :).catch_trial;
                        new_data.ctch(i, :) = temp_catch_trial(trialorder);

%                         new_data.ctch_p1(i, :) = data(mask2, obj.idx.p1);
% 
%                         new_data.ctch_p2(i, :) = data(mask2, obj.idx.p2);
% 
%                         new_data.ctch_corr(i, :) = data(mask2, obj.idx.corr);

%                         temp_cont1 = data(mask, obj.idx.cont1);
%                         new_data.cont1(i, :) = temp_cont1(trialorder);

                        temp_ev2 = data(mask,:).ev2;
                        new_data.ev2(i, :) = temp_ev2(trialorder);

%                         temp_cont2 = data(mask, obj.idx.cont2);
%                         new_data.cont2(i, :) = temp_cont2(trialorder);

                        temp_p1 = data(mask,:).p1;
                        new_data.p1(i, :) = temp_p1(trialorder);

                        temp_p2 = data(mask, :).p2;
                        new_data.p2(i, :) = temp_p2(trialorder);

%                         temp_dist = data(mask,:);
%                         new_data.dist(i, :) = temp_dist(trialorder)./100;

                        temp_rtime = data(mask, :).rtime;
                        new_data.rtime(i, :) = temp_rtime(trialorder);
% 
                        temp_trial = data(mask,:).trial;
                        new_data.trial(i, :) = temp_trial(trialorder);

%                         new_data.catch_trial(i, :) = data(mask2, obj.idx.catch_trial);                        

                        i = i + 1;

                    catch e
                        fprintf(1, '\n There has been an error while treating subject %d \n', i);
                        fprintf(1,'The identifier was:\n%s',e.identifier);
                        fprintf(1,'There was an error! The message was:\n%s',e.message);
                    end
                end

                structlist(isess) = new_data;

            end

            if length(session) > 1
                new_data = obj.merge2sess(structlist);
            end

        end

      

        function new_data = extract_EA(obj, exp_num)
            [data, sub_ids, session, name,nsub] = prepare(obj, exp_num);
            new_data = struct();

            new_data.sess = session;
            new_data.name = name;
            new_data.nsub = nsub;
            new_data.exp_num = exp_num;


            for isess = 1:length(session)
                i = 1;

                for id = 1:length(sub_ids)
                    try
                        sub = sub_ids(id);
                        
                        mask_eli = strcmp(data.phase, 'EA');
                        mask_sub = data.sub_id == sub;
                        mask_catch_trial = data.catch_trial == 0;
%                         mask_ycatch_trial = data(:, obj.idx.catch_trial) == 1;

                        % before exp. 5 op2 has value -1 in ED while after it
                        % takes value 0 (because there were no EE befor                   

                        mask_sess = data.sess ==  session(isess);
                        mask = logical(mask_sub .* mask_sess .* mask_eli .* mask_catch_trial);

                        
%                         mask2 = logical(mask_sub .* mask_sess .* mask_eli .* mask_vs_lot .* mask_ycatch_trial.* mask_vs_lot);

                        trialorder = data(mask,:).trial;

                        if ~issorted(trialorder)
                            [noneed, trialorder] = sort(data(mask,:).trial);
                        else
                            trialorder = 1:length(trialorder);
                        end


                        temp_corr = data(mask, :).corr;
                        new_data.corr(i, :) = temp_corr(trialorder);

                        temp_cho = data(mask, :).cho;
                        new_data.cho(i, :) = temp_cho(trialorder);

                        new_data.cfcho(i, :) = 3 - new_data.cho(i, :);

                        temp_out = data(mask, :).out;
                        new_data.out(i, :) = temp_out(trialorder);

                        temp_ev1 = data(mask, :).ev1;                       

                        new_data.ev1(i, :) = temp_ev1(trialorder);

                        temp_catch_trial = data(mask, :).catch_trial;
                        new_data.ctch(i, :) = temp_catch_trial(trialorder);

%                         new_data.ctch_p1(i, :) = data(mask2, obj.idx.p1);
% 
%                         new_data.ctch_p2(i, :) = data(mask2, obj.idx.p2);
% 
%                         new_data.ctch_corr(i, :) = data(mask2, obj.idx.corr);

%                         temp_cont1 = data(mask, obj.idx.cont1);
%                         new_data.cont1(i, :) = temp_cont1(trialorder);

                        temp_ev2 = data(mask,:).ev2;
                        new_data.ev2(i, :) = temp_ev2(trialorder);

%                         temp_cont2 = data(mask, obj.idx.cont2);
%                         new_data.cont2(i, :) = temp_cont2(trialorder);

                        temp_p1 = data(mask,:).p1;
                        new_data.p1(i, :) = temp_p1(trialorder);

                        temp_p2 = data(mask, :).p2;
                        new_data.p2(i, :) = temp_p2(trialorder);

%                         temp_dist = data(mask,:);
%                         new_data.dist(i, :) = temp_dist(trialorder)./100;

                        temp_rtime = data(mask, :).rtime;
                        new_data.rtime(i, :) = temp_rtime(trialorder);
% 
                        temp_trial = data(mask,:).trial;
                        new_data.trial(i, :) = temp_trial(trialorder);

%                         new_data.catch_trial(i, :) = data(mask2, obj.idx.catch_trial);                        

                        i = i + 1;

                    catch_trial e
                        fprintf(1, '\n There has been an error while treating subject %d \n', i);
                        fprintf(1,'The identifier was:\n%s',e.identifier);
                        fprintf(1,'There was an error! The message was:\n%s',e.message);
                    end
                end


                structlist(isess) = new_data;

            end

            if length(session) > 1
                new_data = obj.merge2sess(structlist);
            end

        end
        function new_data = extract_SA(obj, exp_num)
            [data, sub_ids, session, name,nsub] = prepare(obj, exp_num);
            new_data = struct();

            new_data.sess = session;
            new_data.name = name;
            new_data.nsub = nsub;
            new_data.exp_num = exp_num;


            for isess = 1:length(session)
                i = 1;

                for id = 1:length(sub_ids)
                    try
                        sub = sub_ids(id);

                        mask_eli = data(:, obj.idx.elic) == 0;
                        mask_sub = data(:, obj.idx.sub) == sub;
                        mask_catch_trial = data(:, obj.idx.catch_trial) == 0;
                        mask_vs_amb = data(:, obj.idx.op2) == 2;
                        mask_vs_sym = data(:, obj.idx.op1) == 0;

                        mask_sess = data(:, obj.idx.sess) ==  session(isess);
                        mask = logical(mask_sub .* mask_sess .* mask_eli .* mask_catch_trial .* mask_vs_sym .* mask_vs_amb);

                        trialorder = data(mask, obj.idx.trial);

                        if ~issorted(trialorder)
                            [noneed, trialorder] = sort(data(mask, obj.idx.trial));
                        else
                            trialorder = 1:length(trialorder);
                        end

                        %data = obj.randomize(data, mask, trialorder);


                        temp_corr = data(mask, obj.idx.corr);
                        new_data.corr(i, :) = temp_corr(trialorder);

                        temp_cho = data(mask, obj.idx.cho);
                        new_data.cho(i, :) = temp_cho(trialorder);

                        new_data.cfcho(i, :) = 3 - new_data.cho(i, :);

                        temp_out = data(mask, obj.idx.out);
                        new_data.out(i, :) = temp_out(trialorder);

                        temp_ev1 = data(mask, obj.idx.ev1);
                        new_data.ev1(i, :) = temp_ev1(trialorder);

                        temp_catch_trial = data(mask, obj.idx.catch_trial);
                        new_data.ctch(i, :) = temp_catch_trial(trialorder);

                        temp_cont1 = data(mask, obj.idx.cont1);
                        new_data.cont1(i, :) = temp_cont1(trialorder);

                        temp_ev2 = data(mask, obj.idx.ev2);
                        new_data.ev2(i, :) = temp_ev2(trialorder);

                        temp_cont2 = data(mask, obj.idx.cont2);
                        new_data.cont2(i, :) = temp_cont2(trialorder);

                        temp_p1 = data(mask, obj.idx.p1);
                        new_data.p1(i, :) = temp_p1(trialorder);

                        temp_p2 = data(mask, obj.idx.p2);
                        new_data.p2(i, :) = temp_p2(trialorder);

                        temp_dist = data(mask, obj.idx.dist);
                        new_data.dist(i, :) = temp_dist(trialorder)./100;

                        temp_rtime = data(mask, obj.idx.rtime);
                        new_data.rtime(i, :) = temp_rtime(trialorder);

                        temp_trial = data(mask, obj.idx.trial);
                        new_data.real_trial(i, :) = temp_trial(trialorder);


                        i = i + 1;

                    catch e
                        fprintf(1, '\n There has been an error while treating subject %d \n', i);
                        fprintf(1,'The identifier was:\n%s',e.identifier);
                        fprintf(1,'There was an error! The message was:\n%s',e.message);
                    end
                end

                structlist(isess) = new_data;

            end

            if length(session) > 1
                new_data = obj.merge2sess(structlist);
            end

        end

        function new_data = extract_SP(obj, exp_num)
            
            [data, sub_ids, session, name,nsub] = prepare(obj, exp_num);
            
            new_data = struct();

            new_data.sess = session;
            new_data.name = name;
            new_data.nsub = nsub;
            new_data.exp_num = exp_num;
         
            i = 1;
            for isess = 1:length(session)
                for id = 1:length(sub_ids)
                    try
                        sub = sub_ids(id);
                        
                        mask_eli = strcmp(data.phase, 'SP');
                        mask_sub = data.sub_id == sub;
    
                        mask_catch_trial = data.catch_trial == 0;
    
                        mask_sess = data.sess ==  session(isess);
                        mask = logical(mask_sub .* mask_sess .* mask_eli .* mask_catch_trial);
                        
                        trialorder = data(mask,:).trial;
    
                        if ~issorted(trialorder)
                            [noneed, trialorder] = sort(data(mask,:).trial);
                        else
                            trialorder = 1:length(trialorder);
                        end
    
                        temp_corr = data(mask, :).corr;
                        new_data.corr(i, :) = temp_corr(trialorder);
    
                        temp_cho = data(mask, :).cho;
                        new_data.cho(i, :) = temp_cho(trialorder);
    
                        new_data.cfcho(i, :) = 3 - new_data.cho(i, :);
    
                        temp_out = data(mask, :).out;
                        new_data.out(i, :) = temp_out(trialorder);
    
                        temp_ev1 = data(mask, :).ev1;
    
                        new_data.ev1(i, :) = temp_ev1(trialorder);
    
                        temp_catch_trial = data(mask, :).catch_trial;
                        new_data.ctch(i, :) = temp_catch_trial(trialorder);
    
                  
                        temp_ev2 = data(mask,:).ev2;
                        new_data.ev2(i, :) = temp_ev2(trialorder);
    
    
                        temp_p1 = data(mask,:).p1;
                        new_data.p1(i, :) = temp_p1(trialorder);
    
                        temp_p2 = data(mask, :).p2;
                        new_data.p2(i, :) = temp_p2(trialorder);
    
                        temp_rtime = data(mask, :).rtime;
                        new_data.rtime(i, :) = temp_rtime(trialorder);
                        %
                        temp_trial = data(mask,:).trial;
                        new_data.trial(i, :) = temp_trial(trialorder);
    
    
                        i = i + 1;
                    catch
                        fprintf('There has been an error while treating subject %d \n', i);
                    end
            end
            end
        end
    end

end
