function [data] = get_basic_contacts(data)


contacts = data.contacts; 
NN = data.NNs;

Npop     = data.Npop;
Npop(16) = sum(Npop(16:end));
Npop     = Npop(1:16);
pop_props = [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:16))]/sum(Npop);


%% COMMUNITY-COMMUNITY MATRIX:

C16       = contacts.CM;
Ci        = [C16(:,1),sum(C16(:,2:4),2),sum(C16(:,5:13),2),sum(C16(:,14:16),2)];%sum of the columns
C4        = [Ci(1,:);
            Npop(2:4)'*Ci(2:4,:)/sum(Npop(2:4));
            Npop(5:13)'*Ci(5:13,:)/sum(Npop(5:13));
            Npop(14:16)'*Ci(14:16,:)/sum(Npop(14:16))]; %weighted average of the rows
N4 = [Npop(1), sum(Npop(2:4)), sum(Npop(5:13)), sum(Npop(14:16))];        
Cav      = pop_props*sum(C4,2);
% C4        = contacts.comm*(C4/Cav);
% Cav      = pop_props*sum(C4,2);
contact_props = C4(3,:)/sum(C4(3,:));


%% indices

adInd    = 3;%Adult index
lx       = length(data.obj);%Number of sectors
ln       = length(NN);

NNrel = NN([1:lx,lx+adInd])/sum(NN([1:lx,lx+adInd]));%adult population proportion vector
NNrepvecweighted = zeros(1,ln);
NNrepvecweighted([1:lx,lx+adInd]) = NNrel*contact_props(3);
NNrepvecweighted(lx+[1,2,4]) = contact_props([1,2,4]);
NNrep = repmat(NNrepvecweighted,ln,1);%total population proportion matrix
NNrea = repmat(NN(1:lx)'/sum(NN(1:lx)),lx,1);%workforce population proportion matrix

%% WORKER-WORKER AND COMMUNITY-WORKER MATRICES:

valB          = contacts.B;
valB(lx+1:ln) = 0;
matB          = diag(valB');


valC          = contacts.C;
valC(lx+1:ln) = 0;
matC          = repmat(valC',1,ln).*NNrep;

% move school contacts to students
teacher_contacts = valC(data.EdInd);
frac_infant = data.NNs(lx+1)/sum(data.NNs(lx+[1:2]));
matC(data.EdInd,:) = 0.1 * teacher_contacts .* NNrep(data.EdInd,:);
matC(data.EdInd,lx+1) = matC(data.EdInd,lx+1) + 0.9*frac_infant * teacher_contacts;
matC(data.EdInd,lx+2) = matC(data.EdInd,lx+2) + 0.9*(1-frac_infant) * teacher_contacts;

contacts_between_workers_and_customers = matC .* repmat(NN,1,ln);
contacts_from_customers_to_workers = contacts_between_workers_and_customers' ./ repmat(NN,1,ln);

wnorm = dot(sum(matB+matC+contacts_from_customers_to_workers,2),NN)/sum(NN([1:lx,lx+adInd]));
contacts.workrel = 1;%contacts.workp / wnorm;

matB = contacts.workrel * matB;
av_worker_contacts = dot(sum(matB,2),NN)/sum(NN([1:lx,lx+adInd]));

matC = contacts.workrel * matC;

rel_matC = NNrel' * matC(1:(lx+1),:);
C_distributed = [rel_matC(:,lx+1), rel_matC(:,lx+2), sum(rel_matC(:,[1:lx,lx+3])), rel_matC(:,lx+4)];

contacts_from_customers_to_workers = contacts.workrel * contacts_from_customers_to_workers;
C_back = [sum(contacts_from_customers_to_workers(46,:)), ...
    sum(contacts_from_customers_to_workers(47,:)),...
    sum(NNrel' * contacts_from_customers_to_workers([1:45,48],:)),...
    sum(contacts_from_customers_to_workers(49,:))];

C_back = C_distributed*N4(3) ./ N4;
        
%% get new contact rates
contacts.schoolA1 = C4(1,1) * contacts.schoolA1_frac;
contacts.schoolA2 = C4(2,2) * contacts.schoolA2_frac;

% normalise France values: 18% travel is pt, and pt has 2.5% of contacts,
% which is 0.555
contacts.travelA3 =  contacts.pt/0.18*0.025*Cav;


%% subtract contacts from C4

% school
C4(1,1) = C4(1,1) - contacts.schoolA1;
C4(2,2) = C4(2,2) - contacts.schoolA2;
% travel and work
C4(3,3) = C4(3,3) - contacts.travelA3 * sum(data.NNs(1:45))/sum(data.NNs([1:45,48])) - av_worker_contacts;
% worker to customer
C4(3,:) = C4(3,:) - C_distributed;
% customer to worker
C4(:,3) = C4(:,3) - C_back';

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
C4 = max((1-contacts.hospitality_frac) * C4, 0);

%% save

contacts.C4 = C4;
contacts.contact_props = contact_props;
data.contacts = contacts; 

end