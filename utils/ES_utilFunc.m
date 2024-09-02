function y = ES_utilFunc(betax)

global X Y p_sym p_lot;
% 1-Symbol, 0-lot
y=0;
epsilon=0.00000000001;

betaPre  = betax(1);
midPointPre = betax(2:end);
beta = 600/(1 + exp(-betaPre(1)));
midPoint = 1./(1 + exp(-midPointPre(1:end)));

xUtil1 = [];
xUtil2 = [];
xChoice= [];
for ii=1:size(X,1) %p_sym
    for ij=1:size(X,2) %p_lot
        xUtil1 = [xUtil1; p_lot(ij)]; 
        xUtil2 = [xUtil2; midPoint(ii)]; 
        xChoice= [xChoice; Y(ii,ij)==0];
    end
end

roundsx = length(xChoice);
%find the probability
phat = 1./(1+exp(-beta*(xUtil1-xUtil2)));

phatChosen = zeros(roundsx,1);
for i=1:roundsx
    if(xChoice(i)==1) 
        phatChosen(i) = phat(i);
    else
        phatChosen(i) = 1-phat(i);
    end
end

phatChosen=max(epsilon,phatChosen);
y=y-sum(log(phatChosen));