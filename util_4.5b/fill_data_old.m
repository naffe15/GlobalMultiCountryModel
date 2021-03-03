function yfill = fill_data(y,ya,type,ismooth,iextrap,scale2annual,nostrict)
% yfill = fill_data(y,ya,type)
% 
% y: quarterly series
% ya: annual series
% type = 0 stock
% type = 1 flow
% type = 2 deflator
% ismooth=1: smooth constrained interpolation 
% ismooth=0: non-smooth constrained interpolation
% iextrap=1: extrapolation 
% iextrap=0: no extrapolation
% scale2annual=1: rescale Q data to fit Annual numbers 
% scale2annual=0: rescale A data to fit annualized Q numbers
% nostrict=1: re-scale only if scaling is not too big
% nostrict=0: re-scale indiscriminately with warning message

if nargin==0,
    disp('yfill = fill_data(y,ya,type)')
    return
end

if nargin<4 || isempty(ismooth), ismooth=1; end
if nargin<5 || isempty(iextrap), iextrap=0; end
if nargin<6 || isempty(scale2annual), scale2annual=0; end
if nargin<7 || isempty(nostrict), nostrict=1; end

y=y(:);
ya=ya(:);
yinit = y;
maxa=max(find(~isnan(ya)))*4;
maxq=max(find(~isnan(y)));
if isempty(maxq), maxq=0; end
maxqa=max(maxa,maxq);
y=y(1:maxqa);
ya=ya(1:ceil(maxqa/4));
[nr, ns] = size(y);
[nra, nsa] = size(ya);

if ns>1 || nsa>1,
    error('number of quarterly and annual series must be equal to 1')
end

if nr>nra*4,
    na = ceil(nr/4);
    ya = [ya; nan(na-nra,1)];
    nra=na;
else
    na=nra;
end
if nra*4>nr
    y=[y; nan(nra*4-nr,1)];   %changed by Massimo
    nr = nra*4;
    maxqa=length(y);
end

M0 = eye(nr); % quarterly data are full

switch type
    
    case 0 % stock
        M = eye(nr);
        
    case 1 % flow
        M = ones(nr,nr);
        M=tril(M);
        M=triu(M,-3);
        
    case 2 % deflator
        M = ones(nr,nr);
        M=tril(M);
        M=triu(M,-3);
        M = M./4;
end

inan = find(isnan(y));
if ismember(1,inan) && any(~isnan(y)),
    inan0 = find((1:length(inan))'==inan);
    Mt = transpose(M);
    M(inan0,:) = Mt(inan0,:);
else
    inan0=[];
end

yqa = nan(nr,1);
yqa(4:4:end) = ya;
yqqa = nan(nr,1);
switch type
    case 0 % stock
        yqqa = y;
        
    case 1 % flow
        yqqa = y+lagged(y)+lagged(y,2)+lagged(y,3);
        
    case 2 % deflator
        yqqa = (y+lagged(y)+lagged(y,2)+lagged(y,3))/4;
end

% yscale = mean_nan(yqqa./yqa);
% if ~isnan(yscale),
%    yqa = yqa.*yscale;
% end

yscale = yqqa(4:4:end)./yqa(4:4:end);
if length(find(~isnan(yscale)))>4,
    yscale=interp1(find(~isnan(yscale)), yscale(find(~isnan(yscale))),[1:na]','nearest','extrap');
    if scale2annual && (nostrict==0 || (max(abs(yscale-1))<0.1))
        if max(abs(yscale-1))>=0.1
            disp(['fill_data:: Average scaling:' num2str(mean(abs(yscale-1)))])
        end
        yscale=repmat(yscale',[4 1]);
        yqqa = yqqa./yscale(:);
        y = y./yscale(:);
    else
        if scale2annual
            disp(['WARNING::fill_data:: scaling is not performed since yscale differs from 1 by more than 10%'])
            disp(['fill_data:: MAX scaling:' num2str(max(abs(yscale-1)))])
        end
        yqa(4:4:end) = yqa(4:4:end).*yscale;
    end
end


V = yqa;
for j=1:nr,
    if isnan(y(j)), % missing value to be filled
        M0(j,:)=M(j,:); % re-set observation equation
    end
    if ~isnan(yqqa(j))
        V(j)=yqqa(j);
    end 
end

indx = [1:nr]';
indx0 = indx(~isnan(V));
if any(~isnan(y)) && ismooth==0,
    switch type
        case 0
            yfill = interp1(indx0,V(indx0),indx,'pchip');
        otherwise
            
            Vq = interp1(indx0,V(indx0),indx,'pchip');
            Vqq = y;
            indx1 = indx(isnan(Vqq));
            Vqq(indx1) = Vq(indx1);
            if type,
                Vqq(inan0) = Vq(inan0+3);
            end
            indx2 = indx(~isnan(Vqq));
            yfill = y;
            yfill(indx2) = M0(indx2,indx2)\Vqq(indx2);
    end
    
else % smooth interpolation or no quarterly data available

    ixo=find(~isnan(y));
    istart=(ceil(min(ixo)/4)-1)*4+1;
    is_two_step=0;
    if mod(min(ixo),4)~=1,
        is_two_step=1;
    end
    iend=(ceil(max(ixo)/4))*4;
    if mod(max(ixo),4)~=0,
        is_two_step=1;
    end
    
    switch type
        case 0 %stock
            if iextrap,
            yfill = interp1(indx0,V(indx0),indx,'pchip','extrap');
            else
                yfill=nan(length(indx),1);
                indxa = max(1,indx0(find(~isnan(indx0), 1 ))-3):indx0(find(~isnan(indx0), 1, 'last' ));
                yfill(indxa) = interp1(indx0,V(indx0),indxa,'pchip','extrap');
            end        
        case 1 % flow
            if is_two_step, % incomplete years
                yfill=y;
                yfill(istart:iend) = interp1(indx0-1.5,V(indx0),istart:iend,'pchip','extrap');
                yfill = yfill/4;
                yfill(~isnan(y)) = y(~isnan(y));
                for j=1:nra,
                    ybit = y(1+(j-1)*4:4*j);
                    inox = find(~isnan(ybit));
                    if any(inox) && length(inox)<4,
                        vv=V(4*j)-sum(ybit(inox));
                        inx = (j-1)*4+find(isnan(ybit));
                        yv = sum(yfill( inx ));
                        yscale=vv/yv;
                        yfill(inx)=yfill(inx)*yscale;
                    end
                end
                y = yfill;
                yqqa = [yfill+lagged(yfill)+lagged(yfill,2)+lagged(yfill,3)];
%            V(~isnan(yqqa)) = yqqa(~isnan(yqqa));
            end
            V0=V;
            V=V/4;
            V(~isnan(y)) = y(~isnan(y));
            indx0 = indx(~isnan(V));
     
            if iextrap,
            yfill = interp1(indx0-1.5,V(indx0),indx,'pchip','extrap');
            else
                yfill=nan(length(indx),1);
                indxa = max(1,indx0(find(~isnan(indx0), 1 ))-3):indx0(find(~isnan(indx0), 1, 'last' ));
                yfill(indxa) = interp1(indx0-1.5,V(indx0),indxa,'pchip','extrap');
            end
            %yfill = yfill/4;
            V=V0;
            V(~isnan(yqqa)) = yqqa(~isnan(yqqa));
            yqqa = [yfill+lagged(yfill)+lagged(yfill,2)+lagged(yfill,3)];
            yscale = (V(4:4:end)./yqqa(4:4:end));
            if any(~isnan(yscale)),
               yfill(1:4:end) = yfill(1:4:end).*yscale;
               yfill(2:4:end) = yfill(2:4:end).*yscale;
               yfill(3:4:end) = yfill(3:4:end).*yscale;
               yfill(4:4:end) = yfill(4:4:end).*yscale;
            end
            yfill(~isnan(y)) = y(~isnan(y));
            
        case 2 %deflator

            if is_two_step, % incomplete years
                yfill=y;
                yfill(istart:iend) = interp1(indx0-1.5,V(indx0),istart:iend,'pchip','extrap');
                yfill(~isnan(y)) = y(~isnan(y));
                for j=1:nra,
                    ybit = y(1+(j-1)*4:4*j);
                    inox = find(~isnan(ybit));
                    if any(inox) && length(inox)<4,
                        vv=V(4*j)-sum(ybit(inox))/4;
                        inx = (j-1)*4+find(isnan(ybit));
                        yv = sum(yfill( inx ))/4;
                        yscale=vv/yv;
                        yfill(inx)=yfill(inx)*yscale;
                    end
                end
            y = yfill;
            yqqa = [yfill+lagged(yfill)+lagged(yfill,2)+lagged(yfill,3)]/4;
            end
            V0=V;
            V(~isnan(y)) = y(~isnan(y));
            indx0 = indx(~isnan(V));
            
            if iextrap
            yfill = interp1(indx0-1.5,V(indx0),indx,'pchip','extrap');
            else
                yfill=nan(length(indx),1);
                indxa = max(1,indx0(find(~isnan(indx0), 1 ))-3):indx0(find(~isnan(indx0), 1, 'last' ));
                yfill(indxa) = interp1(indx0-1.5,V(indx0),indxa,'pchip','extrap');
            end
            V=V0;
            V(~isnan(yqqa)) = yqqa(~isnan(yqqa));
            yqqa = [yfill+lagged(yfill)+lagged(yfill,2)+lagged(yfill,3)]/4;
            yscale = (V(4:4:end)./yqqa(4:4:end));
            if any(~isnan(yscale)),
               yfill(1:4:end) = yfill(1:4:end).*yscale;
               yfill(2:4:end) = yfill(2:4:end).*yscale;
               yfill(3:4:end) = yfill(3:4:end).*yscale;
               yfill(4:4:end) = yfill(4:4:end).*yscale;
            end
            yfill(~isnan(y)) = y(~isnan(y));

    end            
end
yfill = [yfill; yinit(maxqa+1:end)];