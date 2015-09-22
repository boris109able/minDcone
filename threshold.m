function[x]=threshold(x,b)
%% thresholding function
for i=1:numel(x)
    if abs(x(i))>b
        x(i) = sign(x(i))*b;
    end 
end
end