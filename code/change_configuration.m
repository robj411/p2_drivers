function new_config = change_configuration(NNvec,lx,i,inext,old_config)
    Xh2w                   = NNvec(1:lx,inext) - NNvec(1:lx,i); %Addition to each wp next intervention step
    Xw2h                   = -Xh2w; 
    Xw2h(Xw2h<0)           = 0;
    Xw2h                   = Xw2h./NNvec(1:lx,i);
    Xw2h(NNvec(1:lx,i)==0) = 0;

    Xh2w(Xh2w<0) = 0;
    Xh2w         = Xh2w/NNvec(lx+3,i);
    

    %Move all infection statuses:
    y0w2h = old_config(1:lx).* Xw2h; % IC%number of people to be put at home (+)
    y0w2h = [-y0w2h; sum(y0w2h,1)];

    y0h2w = old_config(lx+3);
    y0h2w = kron(y0h2w,Xh2w);
    y0h2w = [y0h2w;-sum(y0h2w,1)];
    
    new_config = old_config;
%     new_config([1:lx,lx+3]) = old_config([1:lx,lx+3]) + y0w2h + y0h2w;
end

