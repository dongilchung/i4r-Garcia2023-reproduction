%-------------------------------------------------------------------------
init;
show_current_script_name(mfilename('fullpath'));
%-------------------------------------------------------------------------
selected_exp = [5, 6.1, 6.2, 7.1 ,7.2, 8.1, 8.2, 9.1, 9.2];
save_name = ['data/fit/', 'midpoints_EE_%s_session_%d'];
displayfig = 'on';
force = true;
num = 0;

repl_beta  = [];
repl_midPoint = [];
repl_exp   = [];
repl_subj  = [];
for exp_num = selected_exp
    exp_num
    num = num+ 1;

    sess =  de.get_sess_from_exp_num(exp_num);
    
    data = de.extract_EE(exp_num);
    % ---------------------------------------------------------------------
    % Compute for each subject, the probability of choosing one experienced cue
    %  of choosingsing depending on  cue value
    % --------------------------------------------------------------------
    p_sym = unique(data.p1)';
    nsub = size(data.cho,1);
    
    chose_experience = nan(nsub, length(p_sym), length(p_sym));
    for i = 1:nsub
        for j = 1:length(p_sym)
            count = 0;
            for k = 1:length(p_sym)
                if j ~= k
                    count = count + 1;
                    temp = ...
                        data.cho(...
                            i, logical((data.p2(i, :) == p_sym(k))...
                        .* (data.p1(i, :) == p_sym(j))));
                     chose_experience(i, j, k) = temp == 1;
                end
            end
        end
    end

    for sub = 1:nsub
        sub

        X = zeros(length(p_sym), length(p_sym)-1);
        Y = zeros(length(p_sym), length(p_sym)-1);

        for i = 1:length(p_sym)
            rm = 1:length(p_sym);
            y = reshape(chose_experience(sub, i, :), [], 1);
            y = y(rm~=i);
            x = p_sym;
            x = x(rm~=i);
            Y(i, :) = y;
            X(i, :) = x;
        end

        global X Y p_sym ;
        idx = 1;
        MLE = [];

        for betaPre0 = -15:.5:5
            betaInit=[betaPre0, ones(1, length(p_sym)) .* .5];
            [betaHat,fval,exitflag] = fminsearch(@EE_utilFunc, betaInit, optimset('MaxIter',25000,'MaxFunEvals',25000));
            betaPost  = 1.6*10^6/(1 + exp(-betaHat(1)));
            midPointPost = 1./(1 + exp(-betaHat(2:end)));

            MLE(idx,:) = [exitflag, fval, betaPost midPointPost];

            idx=idx+1;
        end
        raDatabase=[];
        temp=MLE((MLE(:,2)==min(MLE(:,2))),2:end);
        sizeTemp=size(temp,1);
        if(sizeTemp==1)
            % % second 1, original number of solution from the estimation
            raDatabase=[mean(temp,1) 1 NaN];
        else
            temp2=round(temp*10^4)/10^4; % differences smaller than 10^-4 is ignored
            for tempIndex1=1:sizeTemp-1
                for tempIndex2=tempIndex1+1:sizeTemp
                    if sum(temp2(tempIndex1,:)-temp2(tempIndex2,:))~=0 % show 10^-4 difference
                        raDatabase=[mean(temp,1) size(temp,1) NaN]; % show mulitiple local minimum
                        break
                    else
                        raDatabase=[temp2(tempIndex1,:) 1 sizeTemp];
                        % fval rho lamda
                    end
                end
            end
        end
        beta = raDatabase(:,2);
        midPoint = raDatabase(:,3:10);
        repl_beta  = [repl_beta; beta];
        repl_midPoint  = [repl_midPoint; midPoint];
        repl_exp   = [repl_exp; exp_num];
        repl_subj  = [repl_subj; sub];
    end
end
save('data/repl_fit/mh_fit_EE','repl_*')
