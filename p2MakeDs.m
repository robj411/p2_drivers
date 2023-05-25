function [f,data]=p2MakeDs(data,NN,x,hw)

%% COMMUNITY-COMMUNITY MATRIX:

C        = data.CM;
Npop     = data.Npop;
Npop(16) = sum(Npop(16:end));
Npop     = Npop(1:16);

C        = [C(:,1),sum(C(:,2:4),2),sum(C(:,5:13),2),sum(C(:,14:16),2)];%sum of the columns
C        = [C(1,:);
            Npop(2:4)'*C(2:4,:)/sum(Npop(2:4));
            Npop(5:13)'*C(5:13,:)/sum(Npop(5:13));
            Npop(14:16)'*C(14:16,:)/sum(Npop(14:16))];%weighted average of the rows

Cav      = ([Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:16))]/sum(Npop))*sum(C,2);
C        = data.comm*(C/Cav);

%%

adInd    = 3;%Adult index
CworkRow = C(adInd,:);
lx       = length(x);%Number of sectors
ln       = length(NN);

NNrep = repmat(NN'/sum(NN),ln,1);%total population proportion matrix
NNrel = NN([1:lx,lx+adInd])/sum(NN([1:lx,lx+adInd]));%adult population proportion vector
NNrea = repmat(NN(1:lx)'/sum(NN(1:lx)),lx,1);%workforce population proportion matrix

%Make A:
matA                    = zeros(ln,ln);
matA(lx+1:end,lx+1:end) = C;
matA(1:lx,lx+1:end)     = repmat(CworkRow,lx,1);
matA(:,[1:lx,lx+adInd]) = repmat(matA(:,lx+adInd),1,lx+1).*repmat(NNrel',ln,1);

%%

data.EdInd    = 41;%education sector index
data.HospInd  = [32,43,44];%hospitality sector indices

w             = x.^(1/data.alp);
w(data.EdInd) = x(data.EdInd);

if lx==45
    %Education:
    matA(lx+1,lx+1)=    matA(lx+1,lx+1)+    w(data.EdInd)^2*  data.schoolA1;%mixing within age groups only
    matA(lx+2,lx+2)=    matA(lx+2,lx+2)+    w(data.EdInd)^2*  data.schoolA2;
    
    %Hospitality:
    psub=data.NNs(data.HospInd);
    psub=sum(psub.*x(data.HospInd))/sum(psub);%constant from 0-1, weighted measure of how much sectors are open
    matA([1:lx,lx+adInd],:)=    matA([1:lx,lx+adInd],:)+    psub^2*   data.hospA3*    NNrep([1:lx,lx+adInd],:);%mixing between all age groups, including pre-school
    matA(lx+2,:)=               matA(lx+2,:)+               psub^2*   data.hospA2*    NNrep(lx+2,:);
    matA(ln,:)=                 matA(ln,:)+                 psub^2*   data.hospA4*    NNrep(ln,:);

else
    error('Unknown economic configuration!');
    
end

%Transport:
matA(1:lx,1:lx)=    matA(1:lx,1:lx)+    repmat(w',lx,1).*   data.travelA3(1).*  NNrea.*  repmat(1-hw,lx,1).*repmat(1-hw',1,lx);%home-working has a compound effect

%% WORKER-WORKER AND COMMUNITY-WORKER MATRICES:

%Make B and C:
valB          = data.B;
valB          = valB.*(1-hw).*(1-hw);%home-working has a compound effect
valC          = data.C;
valC          = valC.*(1-hw);
valB(lx+1:ln) = 0;
valC(lx+1:ln) = 0;
x(lx+1:ln)    = 0;
w(lx+1:ln)    = 0;
matB          = diag(w.*valB');
matC          = repmat(x.*valC',1,ln).*NNrep;

%%

if ~isfield(data,'wnorm');
    data.wnorm = dot(sum(matB+matC,2),NN)/sum(NN([1:lx,lx+adInd]));
end

D = matA + (data.workp/data.wnorm)*(matB + matC);
f = D;

end