function samples = sample_uniform(variables,CD,country_indices)


samples = zeros(size(variables));
quantile = unifrnd(0,1,1,1);

for i = 1:numel(variables)
    thisone = variables{i};
    values = CD.(thisone);
    values = rmmissing(values(country_indices));
    minval = min(values);
    maxval = max(values);
    samples(i) = unifinv(quantile,minval,maxval);
end

end
