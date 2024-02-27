% moves people aged 20 to 64 to compartments following a change in
% configuration

% NNvec: matrix of configurations
% nSectors: number of sectors
% i: current configuration index
% inext: next configuration index
% old_config: numbers of people in strata 

% new_config: new numbers of people in strata 

function new_config = change_configuration(NNvec,nSectors,i,inext,old_config)
    Xh2w                   = NNvec(1:nSectors,inext) - NNvec(1:nSectors,i); %Addition to each wp next intervention step
    Xw2h                   = -Xh2w; 
    Xw2h(Xw2h<0)           = 0;
    Xw2h                   = Xw2h./NNvec(1:nSectors,i);
    Xw2h(NNvec(1:nSectors,i)==0) = 0;

    Xh2w(Xh2w<0) = 0;
    Xh2w         = Xh2w/NNvec(nSectors+3,i);
    

    %Move all infection statuses:
    y0w2h = old_config(1:nSectors).* Xw2h; %number of people to be put at home (+)
    y0w2h = [-y0w2h; sum(y0w2h,1)];

    y0h2w = old_config(nSectors+3);
    y0h2w = kron(y0h2w,Xh2w);
    y0h2w = [y0h2w;-sum(y0h2w,1)];
    
    new_config = old_config;
%     new_config([1:lx,lx+3]) = old_config([1:lx,lx+3]) + y0w2h + y0h2w;
end

