% take population contact matrix and decompose into items that will be
% impacted differently by configurations
%
% data: struct of general model parameters
%
% contacts: struct of contact parameters

function contacts = get_basic_contacts(data, contacts)


NN = data.NNs;

Npop     = data.Npop;
Npop(16) = sum(Npop(16:end));
Npop     = Npop(1:16);
pop_props = [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:16))]/sum(Npop);


%% COMMUNITY-COMMUNITY MATRIX:

CM_16       = contacts.CM;
CM_164        = [CM_16(:,1),sum(CM_16(:,2:4),2),sum(CM_16(:,5:13),2),sum(CM_16(:,14:16),2)];%sum of the columns
CM_4        = [CM_164(1,:);
            Npop(2:4)'*CM_164(2:4,:)/sum(Npop(2:4));
            Npop(5:13)'*CM_164(5:13,:)/sum(Npop(5:13));
            Npop(14:16)'*CM_164(14:16,:)/sum(Npop(14:16))]; %weighted average of the rows
N4 = [Npop(1), sum(Npop(2:4)), sum(Npop(5:13)), sum(Npop(14:16))];        
CMav      = pop_props*sum(CM_4,2);
contact_props = CM_4(3,:)/sum(CM_4(3,:));


%% indices

adInd    = 3;%Adult index
nSectors       = length(data.obj);%Number of sectors
nStrata       = length(NN);
workage_indices = [1:nSectors,nSectors+adInd];

NNrel = NN(workage_indices)/sum(NN(workage_indices));%adult population proportion vector
NNrepvecweighted = zeros(1,nStrata);
NNrepvecweighted(workage_indices) = NNrel*contact_props(3);
NNrepvecweighted(nSectors+[1,2,4]) = contact_props([1,2,4]);
NNrep = repmat(NNrepvecweighted,nStrata,1);%total population proportion matrix
NNrea = repmat(NN(1:nSectors)'/sum(NN(1:nSectors)),nSectors,1);%workforce population proportion matrix

%% WORKER-WORKER AND COMMUNITY-WORKER MATRICES:

workerworker_contacts          = contacts.B;
workerworker_contacts(nSectors+1:nStrata) = 0;
workerworker_mat          = diag(workerworker_contacts');


customer_to_worker          = contacts.C;
customer_to_worker(nSectors+1:nStrata) = 0;
customer_to_worker_mat          = repmat(customer_to_worker',1,nStrata).*NNrep;

% move school contacts to students
EdInd = data.EdInd;
teacher_contacts = customer_to_worker(EdInd);
frac_infant = NN(nSectors+1)/sum(NN(nSectors+[1:2]));
customer_to_worker_mat(EdInd,:) = 0.1 * teacher_contacts .* NNrep(EdInd,:);
customer_to_worker_mat(EdInd,nSectors+1) = customer_to_worker_mat(EdInd,nSectors+1) + 0.9*frac_infant * teacher_contacts;
customer_to_worker_mat(EdInd,nSectors+2) = customer_to_worker_mat(EdInd,nSectors+2) + 0.9*(1-frac_infant) * teacher_contacts;

% get total contacts between customers and workers
contacts_between_workers_and_customers = customer_to_worker_mat .* repmat(NN,1,nStrata);
% get reciprocal contacts
worker_to_customer_mat = contacts_between_workers_and_customers' ./ repmat(NN,1,nStrata);

% normalise
wnorm = dot(sum(workerworker_mat+customer_to_worker_mat+worker_to_customer_mat,2),NN)/sum(NN(workage_indices));
contacts.workrel = contacts.workp / wnorm;

workerworker_mat = contacts.workrel * workerworker_mat;
av_worker_contacts = dot(sum(workerworker_mat,2),NN)/sum(NN(workage_indices));

customer_to_worker_mat = contacts.workrel * customer_to_worker_mat;

% get marginal contacts by age for workers
rel_mat = NNrel' * customer_to_worker_mat(workage_indices,:);
c_to_w_distributed = [rel_mat(:,nSectors+1), rel_mat(:,nSectors+2), sum(rel_mat(:,workage_indices)), rel_mat(:,nSectors+4)];

% worker_to_customer_mat = contacts.workrel * worker_to_customer_mat;
% C_back = [sum(worker_to_customer_mat(46,:)), ...
%     sum(worker_to_customer_mat(47,:)),...
%     sum(NNrel' * worker_to_customer_mat([1:45,48],:)),...
%     sum(worker_to_customer_mat(49,:))];

c_to_w_back = c_to_w_distributed*N4(3) ./ N4;
        
%% get new contact rates
contacts.schoolA1 = CM_4(1,1) * contacts.schoolA1_frac;
contacts.schoolA2 = CM_4(2,2) * contacts.schoolA2_frac;

% normalise France values: 18% travel is pt, and pt has 2.5% of contacts,
% which is 0.555
contacts.travelA3 =  contacts.pt/0.18*0.025*CMav;


%% subtract contacts from C4

% school
CM_4(1,1) = CM_4(1,1) - contacts.schoolA1;
CM_4(2,2) = CM_4(2,2) - contacts.schoolA2;
% travel and work
CM_4(3,3) = CM_4(3,3) - contacts.travelA3 * sum(NN(1:nSectors))/sum(NN(workage_indices)) - av_worker_contacts;
% customer to worker
CM_4(3,:) = CM_4(3,:) - c_to_w_distributed;
% worker to customer
CM_4(:,3) = CM_4(:,3) - c_to_w_back';

% hospitality
% remaining_contacts = sum(C4,2);

% contacts.hospA2 = contacts.hospitality_frac * remaining_contacts(2);
% contacts.hospA3 = contacts.hospitality_frac * remaining_contacts(3);
% contacts.hospA4 = contacts.hospitality_frac * remaining_contacts(4);
% 
% C4(2,:) = C4(2,:) - contacts.hospA2 * contact_props;
% C4(3,:) = C4(3,:) - contacts.hospA3 * contact_props;
% C4(4,:) = C4(4,:) - contacts.hospA4 * contact_props;

%%!! too many work contacts to infants
CM_4 = max((1-contacts.hospitality_frac) * CM_4, 0);

%% save

contacts.CM_4 = CM_4;
contacts.contact_props = contact_props;

end