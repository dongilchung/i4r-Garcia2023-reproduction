%-------------------------------------------------------------------------
init;
show_current_script_name(mfilename('fullpath'));
%-------------------------------------------------------------------------
selected_exp = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
% save_name = ['data/fit/', 'midpoints_ES_%s_session_%d'];

displayfig = 'on';
force = true;

repl_beta  = [];
repl_midPoint = [];
repl_exp   = [];
repl_subj  = [];
force = 1;
for exp_num = selected_exp
    exp_num

    sess =  de.get_sess_from_exp_num(exp_num);

    data = de.extract_ES(exp_num);

    % ---------------------------------------------------------------------
    % Compute for each experiential symbol of chosing depending on described cue value
    % ---------------------------------------------------------------------

    p_sym = unique(data.p1)';
    p_lot = unique(data.p2)';
    nsub = size(data.cho,1);

    chose_symbol = zeros(nsub, length(p_sym), length(p_lot));
    for i = 1:nsub
        for j = 1:length(p_sym)
            for k = 1:length(p_lot)
                temp = ...
                    data.cho(...
                    i, logical((data.p2(i, :) == p_lot(k))...
                    .* (data.p1(i, :) == p_sym(j))));
                chose_symbol(i, j, k) = temp == 1;

            end
        end
    end

    for sub = 1:nsub
        sub
        X = zeros(length(p_sym), length(p_lot));
        Y = zeros(length(p_sym), length(p_lot));

        for i = 1:length(p_sym)
            Y(i, :) = reshape(chose_symbol(sub, i, :), [], 1);
            X(i, :) = p_lot;
        end

        global X Y p_sym p_lot;
        idx = 1;
        MLE = [];

        for betaPre0 = -5:.5:5
            % for midPoint1Pre0 = -5:.5:5
            %     for midPoint2Pre0 = -5:.5:5
            %         for midPoint3Pre0 = -5:.5:5
            %             for midPoint4Pre0 = -5:.5:5
            %                 for midPoint5Pre0 = -5:.5:5
            %                     for midPoint6Pre0 = -5:.5:5
            %                         for midPoint7Pre0 = -5:.5:5
            %                             for midPoint8Pre0 = -5:.5:5
            % betaInit=[betaPre0, midPoint1Pre0,...
            %     midPoint2Pre0,midPoint3Pre0,midPoint4Pre0,midPoint5Pre0,midPoint6Pre0,...
            %     midPoint7Pre0,midPoint8Pre0];
            betaInit=[betaPre0, ones(1, length(p_sym)) .* .5];
            [betaHat,fval,exitflag] = fminsearch(@ES_utilFunc, betaInit, optimset('MaxIter',25000,'MaxFunEvals',25000));
            betaPost  = 600/(1 + exp(-betaHat(1)));
            midPointPost = 1./(1 + exp(-betaHat(2:end)));

            MLE(idx,:) = [exitflag, fval, betaPost midPointPost];

            idx=idx+1;
            %                         end
            %                     end
            %                 end
            %             end
            %         end
            %     end
            % end
            % end
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
        tmp1 = raDatabase(:,3:2+length(p_sym));
        tmp2 = zeros(1,8);
        tmp2(1,1:length(tmp1)) = tmp1;
        midPoint = tmp2;
        repl_beta  = [repl_beta; beta];
        repl_midPoint  = [repl_midPoint; midPoint];
        repl_exp   = [repl_exp; exp_num];
        repl_subj  = [repl_subj; sub];
    end
end
save('data/repl_fit/mh_fit_ES','repl_*')

