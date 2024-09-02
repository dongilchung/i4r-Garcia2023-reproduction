function score = heuristic(data)

for sub = 1:size(data.cho,1)
    
    for t = 1:size(data.cho,2)
        
        
        if data.p2(sub,t) >= .5
            prediction = 2;
        else
            prediction = 1;
        end
        
        score(sub, t) = prediction;
        
    end
end
end