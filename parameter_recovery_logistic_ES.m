%-------------------------------------------------------------------------
init;
show_current_script_name(mfilename('fullpath'));
%-------------------------------------------------------------------------
selected_exp = [1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2];
save_name = ['data/pr/', 'midpoints_ES_%s_session_%d'];

displayfig = 'on';
force = true;

for exp_num = selected_exp

    %---------------------------------------------------------------------%
    % get data parameters                                                           %
    % --------------------------------------------------------------------%
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

    param = load(...
        sprintf('data/fit/midpoints_%s_%s_session_%d',...
        'ES', name, sess));

    pre.midpoints = param.midpoints;
    pre.beta1  = param.beta1;

    fprintf('Fitting exp. %s \n', num2str(exp_num));
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
                xUtil1 = p_lot(k);
                xUtil2 = pre.midpoints(i,j);
                phat = 1./(1+exp(-pre.beta1(i)*(xUtil2-xUtil1)));
                if phat>rand
                    tempCh = 1;
                else
                    tempCh = 0;
                end
                chose_symbol(i, j, k) = tempCh;
            end
        end
    end

    midpoints = nan(nsub, length(p_sym));
    params = nan(nsub, length(p_sym)+1);
    beta1 = nan(nsub, 1);
    nll = nan(nsub, 1);

    for sub = 1:nsub

        X = zeros(length(p_sym), length(p_lot));
        Y = zeros(length(p_sym), length(p_lot));

        for i = 1:length(p_sym)
            Y(i, :) = reshape(chose_symbol(sub, i, :), [], 1);
            X(i, :) = p_lot;
        end

        try
            if force
                error('fitting');
            end
            param = load(...
                sprintf(save_name,...
                data.name, sess ...
                ));
            beta1 = param.beta1;
            midpoints = param.midpoints;
            nll = param.nll;
            tosave = false;
        catch
            tosave = true;
            options = optimset(...
                'Algorithm',...
                'interior-point',...
                'Display', 'off',...
                'MaxIter', 10000,...
                'MaxFunEval', 10000);

            [params(sub, :), nll(sub)] = fmincon(...
                @(x) mle(x, X, Y),...
                [1, ones(1, length(p_sym)) .* .5],...
                [], [], [], [],...
                [0.01, zeros(1, length(p_sym))],...
                [inf, ones(1, length(p_sym))],...
                [],...
                options...
                );

        end

    end

    midpoints = params(:, 2:length(p_sym)+1);
    beta1 = params(:, 1);
    if tosave
        param.midpoints = midpoints;
        param.beta1 = beta1;
        param.nll = nll;

        save(sprintf(save_name, data.name, sess),...
            '-struct', 'param');
    end

end



function nll = mle(params, X, Y)
options = optimset('Display','off');
temp = params(1);
midpoints = params(2:end);
ll = 0;
for i = 1:size(Y, 1)
    yhat = logfun(X(i,:), midpoints(i), temp);
    ll = ll + (1/numel(yhat)) * sum(log(yhat) .* Y(i,:) + log(1-yhat).*(1-Y(i,:)));
end
nll = -ll;
end

function p = logfun(x, midpoint, temp)
p = 1./(1+exp(temp.*(x-midpoint)));
end

