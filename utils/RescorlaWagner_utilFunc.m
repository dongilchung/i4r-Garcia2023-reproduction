function y = RescorlaWagner_utilFunc(betax)

global xData sData;
x=xData;
s=sData;

y=0;
epsilon=0.00000000001;

alphaPre = betax(1);
betaPre  = betax(2);

alpha = 1/(1 + exp(-alphaPre(1)));
beta = 50/(1 + exp(-betaPre(1)));

xChoice   = x.cho(s,:);
roundsx   = size(xChoice,2);
xCon= x.con(s,:);
nCon= length(unique(xCon));
xOut= x.out(s,:);

xOutCf= x.cfout(s,:); % unchosen outcome

Qset1  = zeros(nCon, roundsx);
Qset2  = zeros(nCon, roundsx);

Qset1(:,1) = 0.5;
Qset2(:,1) = 0.5;
Qhist1 = [];
Qhist2 = [];
for trial = 1:roundsx-1
    Qset1(:,trial+1) = Qset1(:,trial);
    Qset2(:,trial+1) = Qset2(:,trial);
    tempQ1 = Qset1(xCon(trial),trial);
    tempQ2 = Qset2(xCon(trial),trial);
    Qhist1 = [Qhist1; tempQ1];
    Qhist2 = [Qhist2; tempQ2];

    % phat = 1./(1+exp(-beta*(tempQ1-tempQ2)));
    % if phat>rand
    %     tempCh = 1;
    % else
    %     tempCh = 2;
    % end

    tempCh   = xChoice(trial);
    chosenR  = xOut(trial);
    unchosenR= xOutCf(trial);

    if x.fit_cf
        if tempCh==1
            tempPEc = chosenR - tempQ1;
            tempPEu = unchosenR - tempQ2;
            tempQ1 = tempQ1 + alpha.*tempPEc;
            tempQ2 = tempQ2 + alpha.*tempPEu;
        else
            tempPEc = chosenR - tempQ2;
            tempPEu = unchosenR - tempQ1;
            tempQ1 = tempQ1 + alpha.*tempPEu;
            tempQ2 = tempQ2 + alpha.*tempPEc;
        end
    else
        if tempCh==1
            tempPEc = chosenR - tempQ1;
            tempPEu = 0;
            tempQ1 = tempQ1 + alpha.*tempPEc;
            tempQ2 = tempQ2 + alpha.*tempPEu;
        else
            tempPEc = chosenR - tempQ2;
            tempPEu = 0;
            tempQ1 = tempQ1 + alpha.*tempPEu;
            tempQ2 = tempQ2 + alpha.*tempPEc;
        end
    end
    Qset1(xCon(trial),trial+1) = tempQ1;
    Qset2(xCon(trial),trial+1) = tempQ2;
    if trial==roundsx-1
        Qhist1 = [Qhist1; tempQ1];
        Qhist2 = [Qhist2; tempQ2];
    end
end

%find the probability
phat = 1./(1+exp(-beta*(Qhist1-Qhist2)));

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