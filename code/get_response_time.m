function response_time = get_response_time(data, dis, hosp_trigger)
    
    Npop = data.Npop;
    ihr = dis.ihr';
    Npop(length(ihr)) = sum(Npop(length(ihr):length(Npop)));
    Npop = Npop(1:length(ihr));
    mean_ihr = sum(Npop .* ihr) / sum(Npop);
    
    R0 = dis.R0; 
    nsamples = 100;
    generationsam = zeros(1,nsamples);
    for i = 1:nsamples
        generationsam(i) = branchingProcess(R0, mean_ihr, hosp_trigger);
    end
    generations = mean(generationsam(generationsam>0));
    
    generation_time = log(R0) / (log(2) / dis.Td);
    response_time = generations * generation_time + dis.Tlat + dis.Tsh;


end


function generations = branchingProcess(R0, mean_ihr, hosp_trigger)
    % Initialize variables
    finalSize = [];
    hospitalisations = 0;
    generations = 1;
    finalSize(generations) = 5;
    while hospitalisations < hosp_trigger & finalSize(generations) > 0
        generations = generations + 1;
        % In subsequent generations, each infected individual generates R0 offspring
        offspring = poissrnd(R0 * finalSize(generations - 1));
        finalSize(generations) = offspring;
        hospitalisations = binornd(finalSize(generations), mean_ihr);
    end
end
