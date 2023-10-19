function f = p2MakeDs(data,NN,x,hw)


contacts = data.contacts; 

w             = x;
w(data.EdInd) = x(data.EdInd);


%% COMMUNITY-COMMUNITY MATRIX:
C4 = contacts.C4;
contact_props = contacts.contact_props;

%% add contacts to C4

% start with hospitality
psub = data.NNs(data.HospInd);
psub = sum(psub.*x(data.HospInd))/sum(psub);%constant from 0-1, weighted measure of how much sectors are open
C4 = C4 + psub^2 * contacts.hospitality_frac / (1-contacts.hospitality_frac) * C4;

% school
C4(1,1) = C4(1,1) + w(data.EdInd).^2 * contacts.schoolA1;
C4(2,2) = C4(2,2) + w(data.EdInd).^2 * contacts.schoolA2;

% hospitality
% C4(2,:) = C4(2,:) + psub^2 * contacts.hospA2 * contact_props;
% C4(3,:) = C4(3,:) + psub^2 * contacts.hospA3 * contact_props;
% C4(4,:) = C4(4,:) + psub^2 * contacts.hospA4 * contact_props;

%%

adInd    = 3;%Adult index
CworkRow = C4(adInd,:);
lx       = length(x);%Number of sectors
ln       = length(NN);

NNrel = NN([1:lx,lx+adInd])/sum(NN([1:lx,lx+adInd]));%adult population proportion vector
NNrepvecweighted = zeros(1,ln);
NNrepvecweighted([1:lx,lx+adInd]) = NNrel*contact_props(3);
NNrepvecweighted(lx+[1,2,4]) = contact_props([1,2,4]);
NNrep = repmat(NNrepvecweighted,ln,1);%total population proportion matrix
% NNrep = repmat(NN'/sum(NN),ln,1);%total population proportion matrix
NNrea = repmat(NN(1:lx)'/sum(NN(1:lx)),lx,1);%workforce population proportion matrix

%Make A:
matA                    = zeros(ln,ln);
matA(lx+1:end,lx+1:end) = C4;
matA(1:lx,lx+1:end)     = repmat(CworkRow,lx,1);
matA(:,[1:lx,lx+adInd]) = repmat(matA(:,lx+adInd),1,lx+1).*repmat(NNrel',ln,1);

%Transport:
matA(1:lx,1:lx)=    matA(1:lx,1:lx)+    repmat(w',lx,1).*   contacts.travelA3(1).*  NNrea.*  repmat(1-hw,lx,1).*repmat(1-hw',1,lx);%home-working has a compound effect
mat = repmat(w',lx,1).*   contacts.travelA3(1).*  NNrea.*  repmat(1-hw,lx,1).*repmat(1-hw',1,lx);

%% WORKER-WORKER AND COMMUNITY-WORKER MATRICES:

%Make B and C:
valB          = contacts.B;
valB          = valB.*(1-hw).*(1-hw);%home-working has a compound effect
valC          = contacts.C;
valC          = valC.*(1-hw);
valB(lx+1:ln) = 0;
valC(lx+1:ln) = 0;
x(lx+1:ln)    = 0;
w(lx+1:ln)    = 0;
matB          = diag(w.*valB');
matC          = repmat(x.*valC',1,ln).*NNrep;
% move school contacts to students
teacher_contacts = x(data.EdInd).*valC(data.EdInd);
matC(data.EdInd,:) = 0.1 * teacher_contacts .* NNrep(data.EdInd,:);
frac_infant = data.NNs(lx+1)/sum(data.NNs(lx+[1:2]));
matC(data.EdInd,lx+1) = matC(data.EdInd,lx+1) + 0.9*frac_infant * teacher_contacts;
matC(data.EdInd,lx+2) = matC(data.EdInd,lx+2) + 0.9*(1-frac_infant) * teacher_contacts;

contacts_between_workers_and_customers = matC .* repmat(NN,1,ln);
Cback = contacts_between_workers_and_customers' ./ repmat(NN,1,ln);

D = matA + contacts.workrel*(matB + matC + Cback);

f = D;

end