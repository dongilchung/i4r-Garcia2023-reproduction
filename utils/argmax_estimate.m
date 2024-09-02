function score = argmax_estimate(data, symp, values)
for sub = 1:size(data.cho,1)
    
    for t = 1:size(data.cho,2)
        
        
        if data.p2(sub,t) >= values(sub, symp==data.p1(sub,t))
            prediction = 2;
        else
            prediction = 1;
        end
        
        score(sub, t) = prediction;
        
    end
end
end
