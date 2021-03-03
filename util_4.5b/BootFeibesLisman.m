function bflOutput=BootFeibesLisman(Y1,d,ToD,OoD);

% It implements Boot Feibes and Lisman Procedure 


% 	--------INPUT------------

% 	Y= low frequencY time series
% 	ToD: type of temporal aggregation 
% 	ToD=1,2,3,4: sum (flow), average (flow), last (stock) and First (stock)
% 	OoD: Order of Disaggregation
% 	OoD=12,4,2,6,2,3: Y2M, Y2Q, Y2HY, HY2M,HY2Q, Q2M
% 	d= order of differences


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%OUTPUT INIZIALIZATION%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bflOutput.y = [];
bflOutput.d  = [];
bflOutput.B       = [];
bflOutput.H       = [];
bflOutput.Y       = [];
bflOutput.OoD      = [];
bflOutput.ToD       = [];

% % % % % % % % % % % % % % % % % % % % % % % % % i=1;

[LfD,NfD]=size(Y1);

switch ToD
    case 1   
      c=ones(1,OoD);
    case 2
    c=ones(1,OoD)./OoD;
    case 3
     c=zeros(1,OoD);
     c(OoD)=1;
    case 4
     c=zeros(1,OoD);
     c(1)=1;
    otherwise
     error ('  Type of disaggregation in not defined ');
end

B=kron(eye(LfD),c);

y=NaN*ones(OoD*LfD,NfD);
for ii=1:NfD
    Y=Y1(:,ii);
    qq=isnan(Y);
    
    [aa,bb]=find(qq==0);
    K=[];
    for jj=1:length(aa)
        k=(OoD*(aa(jj)-1)+1:OoD*aa(jj));
        K=[K;k];
    end
    Y=Y(~isnan(Y));
    
    
    % Number of HF points
    
    n=OoD*LfD; 
    
    % Aggregation matrix 

if d==0
    Temp_D=eye(n);
else
Temp_D=zeros(n-d,n);


v=zeros(d,d+1);
v(:,1:2)=1;

for i=2:d;
    for j=2:i+1
        v(i,j)=v(i-1,j-1)+v(i-1,j);
    end
end

a=zeros(1,d+1);

    for i=1:d+1
        a(i)=(-1)^(i+1);
    end
v=v(end,1:d+1);
v=v.*a;

    for i=1:n-d
        Temp_D(i,i:i+d)=v;
    end
Temp_D=(-1)^d*Temp_D;
end

    
    % Matrix of differentiation
    
 
    S=2*Temp_D'*Temp_D;
    F=zeros(LfD,LfD);
    
    % Filterinng matrix
    
    H=[S -B'; B F];
    Z=[zeros(n,NfD);Y];
    
    % Estimates  
    
    Temp_Y=H\Z;
    y1=Temp_Y(1:n,:);
    % y(end-OoD*sY+1:end,i)=y1;
    y(K,ii)=0;
    y(~isnan(y),ii)=y1;
    
end

bflOutput.data='Estimated Time Series' ;   
bflOutput.y = y;
bflOutput.Amatrix='Aggregation Matrix' ;   
bflOutput.B       = B;
bflOutput.Fmatrix='Filtering  Matrix' ;   
bflOutput.H       = H;
bflOutput.strDoD='Degree of differencing' ;   
bflOutput.d       = d;
bflOutput.LFTS='Aggregated Time Series' ;   
bflOutput.Y       = Y;

bflOutput.OoD      = OoD;
bflOutput.ToD       = ToD;

% end