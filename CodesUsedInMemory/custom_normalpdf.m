function out = custom_normalpdf(x,mu,sigmas)
out = zeros(length(sigmas),length(x));
row_index = 1;
    for sigma = sigmas
        temp = normpdf(x,mu,sigma);
        
        out(row_index,:) = temp;
        
        row_index = row_index + 1;
        
    end
end