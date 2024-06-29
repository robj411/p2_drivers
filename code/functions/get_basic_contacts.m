% take population contact matrix and decompose into items that will be
% impacted differently by configurations
%
% data: struct of general model parameters
%
% contacts: struct of contact parameters

function contacts = get_basic_contacts(data, contacts)


NN = data.NNs;
CM_16       = contacts.CM;
s16 = size(CM_16,1);

Npop     = data.Npop;
Npop(s16) = sum(Npop(s16:end));
Npop     = Npop(1:s16);
ageindex = data.ageindex;
ageindex{4} = min(ageindex{4}):s16;
Npop4 = data.Npop4;
pop_props = Npop4/sum(Npop);


%% COMMUNITY-COMMUNITY MATRIX:

CM_164    = cell2mat(arrayfun(@(x) sum(CM_16(:,x{1}),2), ageindex, 'UniformOutput', false)); %sum of the columns
CM_4      = cell2mat(arrayfun(@(x) (Npop(x{1})'*CM_164(x{1},:)/sum(Npop(x{1})))', ageindex, 'UniformOutput', false))';        
CMav      = pop_props*sum(CM_4,2);
contact_props = CM_4(3,:)/sum(CM_4(3,:));
workage_total = sum(CM_4(3,:));

%% indices

adInd    = data.adInd; %Adult index
nSectors       = length(data.obj); %Number of sectors
nStrata       = length(NN);
workage_indices = [1:nSectors,nSectors+adInd];

NNrel = NN(workage_indices)/sum(NN(workage_indices)); %adult population proportion vector

%% WORKER-WORKER AND COMMUNITY-WORKER MATRICES:

sectoragedist = zeros(nSectors,nStrata);
sectorcontactfracs = contacts.sectorcontactfracs;
% correct for age dist
over65frac = Npop4(4)/sum(Npop4);
ukover65frac = .25;
sectorcontactfracs.x65plus = sectorcontactfracs.x65plus * over65frac / ukover65frac;
newtotal = sectorcontactfracs.workingage+sectorcontactfracs.x65plus+sectorcontactfracs.under18;
sectorcontactfracs.under18 = sectorcontactfracs.under18 ./ newtotal;
sectorcontactfracs.workingage = sectorcontactfracs.workingage ./ newtotal;
sectoragedist(:,nSectors+[1,2]) = sectorcontactfracs.under18 * pop_props(1:2)/sum(pop_props(1:2));
sectoragedist(:,nSectors+4) = sectorcontactfracs.x65plus;
sectoragedist(:,workage_indices) = sectorcontactfracs.workingage * NNrel';

community_to_worker          = contacts.sectorcontacts;
community_to_worker_mat          = repmat(community_to_worker,1,nStrata).*sectoragedist;
community_to_worker_mat(nSectors+1:nStrata,:) = 0;


% normalise
worker_total = dot(sum(community_to_worker_mat,2),NN)/sum(NN(workage_indices));
target_work_contacts = contacts.work_frac*workage_total;
contacts.work_scalar = target_work_contacts / worker_total;

community_to_worker_mat = contacts.work_scalar * community_to_worker_mat;
contacts.community_to_worker_mat = community_to_worker_mat;

% get marginal contacts by age for workers
rel_mat = NNrel' * community_to_worker_mat(workage_indices,:);
c_to_w_distributed = [rel_mat(:,nSectors+1), rel_mat(:,nSectors+2), sum(rel_mat(:,workage_indices)), rel_mat(:,nSectors+4)];

c_to_w_back = c_to_w_distributed*Npop4(3) ./ Npop4;
        
%% get new contact rates
contacts.school1 = CM_4(1,1) * contacts.school1_frac;
contacts.school2 = CM_4(2,2) * contacts.school2_frac;


%% subtract contacts from C4

% school
CM_4(1,1) = CM_4(1,1) - contacts.school1;
CM_4(2,2) = CM_4(2,2) - contacts.school2;
% customer to worker
CM_4(3,:) = CM_4(3,:) - c_to_w_distributed;
% worker to customer
CM_4([1,2,4],3) = CM_4([1,2,4],3) - c_to_w_back([1,2,4])';

% hospitality

%%!! too many work contacts to infants
hospitality_age = contacts.hospitality_age;
hospitality_age = [hospitality_age(1,:); hospitality_age];
total_contacts = sum(CM_4,2);
contacts.hospitality_contacts = repmat(total_contacts .* contacts.hospitality_frac,1,4) .* hospitality_age;
CM_4 = max(CM_4 - contacts.hospitality_contacts, 0);

%% save

contacts = rmfield(contacts,'sectorcontactfracs');
contacts.CM_4 = CM_4;
contacts.contact_props = contact_props;

end



