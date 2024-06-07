% simulates a random response time based on country and pathogen
% characteristics, taking the average of many simulations

% dis: struct of pathogen parameters
% data: struct of general model parameters
% hosp_trigger: number of hospital cases at which the country responds

% response_time: time  at which the country responds
% (international alert time)

function response_time = get_response_time(data, dis, hosp_trigger)
    
    Npop = data.Npop;
    ihr = dis.ihr';
    Npop(length(ihr)) = sum(Npop(length(ihr):length(Npop)));
    Npop = Npop(1:length(ihr));
    mean_ihr = sum(Npop .* ihr) / sum(Npop);
    
    R0 = dis.R0; 
    nsamples = 100;
    generationsam = zeros(1,nsamples);
    hospsam = zeros(1,nsamples);
    for i = 1:nsamples
        branch_out = branching_process(R0, mean_ihr, hosp_trigger);
        generationsam(i) = branch_out(1);
        hospsam(i) = branch_out(2);
    end
    generations = mean(generationsam(hospsam>=hosp_trigger));
    
    response_time = generations * dis.generation_time;


end


% simulates one branching process

% R0: basic reproduction number 
% mean_ihr: average ratio of hospitalisation
% hosp_trigger: number of hospital cases at which a response is made

% branch_out: vector with the number of generations and the number of
% hospitalisations in the simulation

function branch_out = branching_process(R0, mean_ihr, hosp_trigger)
    % Initialize variables
    n_infections = [];
    hospitalisations = 0;
    generations = 1;
    n_infections(generations) = 5;
    while hospitalisations < hosp_trigger && n_infections(generations) > 0
        generations = generations + 1;
        % In subsequent generations, each infected individual generates R0 offspring
        n_infections(generations) = poissrnd(R0 * n_infections(generations - 1));
        hospitalisations = hospitalisations + binornd(n_infections(generations), mean_ihr);
    end
    branch_out = [generations, hospitalisations];
end



