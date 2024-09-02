%-------------------------------------------------------------------------
init;
%-------------------------------------------------------------------------

%% A5. Correct choice rate in the LE phase: Descriptive statistics of learning performance across Experiments 1-4


% ------------------------------------------------------------------------
% customizable parameters
% ------------------------------------------------------------------------
selected_exp = [1, 2, 3, 4];
modality = 'LE';
color = blue;
displayfig = 'on';
filename = 'Fig2A';

% ------------------------------------------------------------------------
% fixed parameters
% ------------------------------------------------------------------------
T_exp = table();
T_con = table();


sub_count = 0;
num = 0;
stat_A1 = [];
for exp_num = selected_exp
    clear dd
    num = num + 1;

    data = de.extract_LE(exp_num);
    data_ed = de.extract_ES(exp_num);

    if exp_num == 4
        data.con(data.con == 2) = 4;
    end
    ncon = length(unique(data.con));

    dd = NaN(ncon, data.nsub);
    cons = flip(unique(data.con));

    for i = 1:ncon
        for sub = 1:data.nsub

            dd(i, sub) = mean(...
                data.corr(sub, data.con(sub,:)==cons(i)));
            %if ismember(cons(i), [1, 4])
            complete = ismember(exp_num, [3, 4]);
            block = ismember(exp_num, [2, 3, 4]);
            less_cues = exp_num == 4;


            T3 = table(...
                sub+sub_count, exp_num,  complete,  block, less_cues, dd(i, sub), cons(i), ...
                'variablenames',...
                {'subject', 'exp_num', 'complete', 'block','less_cues', 'score', 'cond'}...
                );

            T_con = [T_con; T3];
            % end
        end
    end

    for sub = 1:data.nsub
        s1 = mean(data.corr(sub, :));
        s2 = mean(data_ed.corr(sub, :));

        complete = int32(ismember(exp_num, [3, 4]));
        block = int32(ismember(exp_num, [2, 3, 4]));
        less_cues = int32(exp_num == 4);

        T1 = table(...
            sub+sub_count, exp_num, complete,  block, less_cues, s1, {'LE'}, ...
            'variablenames',...
            {'subject', 'exp_num', 'complete', 'block','less_cues', 'score', 'modality'}...
            );

        T_exp = [T_exp; T1];
    end
    addpath(genpath('./repl_utils'))
    idx = 1;
    if exp_num~=4
        num_cond = [1:4];
    else
        num_cond = [1,4];
    end
end

% ------------------------------------------------------------------------

stat_A5 = [];
for exp_num = selected_exp
    y = T_con.score;
    x = T_con.cond;
    z = T_con.exp_num;
    s = T_con.subject;

    tmp_x = x(z==exp_num);
    tmp_y = y(z==exp_num);
    tmp_s = s(z==exp_num);

    stat_A5 = [stat_A5; exp_num mean(tmp_y) std(tmp_y)./sqrt(length(tmp_y))];
end
