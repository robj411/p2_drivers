% construct contact matrices from components for configurations
%
% data: struct of general model parameters
% NN: population by stratum
% x: economic configuration
% hw: proportion working from home by stratum
%
% D: contact matrix

function contact_matrix = p2MakeDs(data,NN,x,hw)

%% variables to use

contacts = data.contacts; 
CM_4 = contacts.C4;
contact_props = contacts.contact_props;
hospitality_frac = contacts.hospitality_frac;

w             = x;
edInd = data.EdInd;
adInd    = data.adInd; %Adult index
HospInd = data.HospInd;

workRow = CM_4(adInd,:);
nSectors       = length(x);%Number of sectors
nStrata       = length(NN);
workage_indices = [1:nSectors,nSectors+adInd];

NNrel = NN(workage_indices)/sum(NN(workage_indices));%adult population proportion vector
NNrepvecweighted = zeros(1,nStrata);
NNrepvecweighted(workage_indices) = NNrel*contact_props(3);
NNrepvecweighted(nSectors+[1,2,4]) = contact_props([1,2,4]);
NNrep = repmat(NNrepvecweighted,nStrata,1);%total population proportion matrix
NNrea = repmat(NN(1:nSectors)'/sum(NN(1:nSectors)),nSectors,1);%workforce population proportion matrix

%% add school and hospitality contacts to CM_4

% start with hospitality
hospitality_sectors = NN(HospInd);
% get weighted average
hospitality_sectors = sum(hospitality_sectors.*x(HospInd))/sum(hospitality_sectors); % constant from 0 to 1, weighted measure of how much sectors are open
CM_4 = CM_4 + hospitality_sectors^2 * hospitality_frac / (1-hospitality_frac) * CM_4;

% school
CM_4(1,1) = CM_4(1,1) + w(edInd).^2 * contacts.schoolA1;
CM_4(2,2) = CM_4(2,2) + w(edInd).^2 * contacts.schoolA2;

%% Make community matrix
community_mat                    = zeros(nStrata,nStrata);
community_mat(nSectors+1:end,nSectors+1:end) = CM_4;
community_mat(1:nSectors,nSectors+1:end)     = repmat(workRow,nSectors,1);
community_mat(:,workage_indices) = repmat(community_mat(:,nSectors+adInd),1,nSectors+1).*repmat(NNrel',nStrata,1);

%Transport:
community_mat(1:nSectors,1:nSectors) = community_mat(1:nSectors,1:nSectors) + ...
    repmat(w',nSectors,1).*repmat(w,1,nSectors).*contacts.travelA3(1).*NNrea.*repmat(1-hw,nSectors,1).*repmat(1-hw',1,nSectors); % home-working has a compound effect
% mat = repmat(w',nSectors,1).*   contacts.travelA3(1).*  NNrea.*  repmat(1-hw,nSectors,1).*repmat(1-hw',1,nSectors);

%% WORKER-WORKER AND COMMUNITY-WORKER MATRICES

x(nSectors+1:nStrata)    = 0;
w(nSectors+1:nStrata)    = 0;

% workers
workerworker_val          = contacts.B;
workerworker_val          = workerworker_val.*(1-hw).*(1-hw);%home-working has a compound effect
workerworker_val(nSectors+1:nStrata) = 0;
workerworker_mat          = diag(w.^2.*workerworker_val');

% customer to worker
customertoworker_val          = contacts.C;
customertoworker_val          = customertoworker_val.*(1-hw).*(1-hw);
customertoworker_val(nSectors+1:nStrata) = 0;
customertoworker_mat          = repmat(x.*w.*customertoworker_val',1,nStrata).*NNrep;

% move school contacts to students
frac_infant = NN(nSectors+1)/sum(NN(nSectors+[1:2]));
teacher_contacts = x(edInd).*customertoworker_val(edInd);
customertoworker_mat(edInd,:) = 0.1 * teacher_contacts .* NNrep(edInd,:);
customertoworker_mat(edInd,nSectors+1) = customertoworker_mat(edInd,nSectors+1) + 0.9*frac_infant * teacher_contacts;
customertoworker_mat(edInd,nSectors+2) = customertoworker_mat(edInd,nSectors+2) + 0.9*(1-frac_infant) * teacher_contacts;

% get return contacts
contacts_between_workers_and_customers = customertoworker_mat .* repmat(NN,1,nStrata);
worker_back = contacts_between_workers_and_customers' ./ repmat(NN,1,nStrata);

%% add all together
contact_matrix = community_mat + contacts.workrel*(workerworker_mat + customertoworker_mat + worker_back);

end