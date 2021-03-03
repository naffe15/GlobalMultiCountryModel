function MAP0=CreateColorMap(nbcmpts)

MAP0=ones(nbcmpts,3);

% MAP1=MAP0;


if nbcmpts == 1
    MAP0(1,:)=[1 0 0];
elseif    nbcmpts == 2
    
    MAP0(1,:)=[0 0 1];
    MAP0(2,:)=[1 0 0];
elseif nbcmpts == 3
    
    MAP0(1,:)=[0 0 1];
    MAP0(2,:)=[0 1 0];
    MAP0(3,:)=[1 0 0];
elseif nbcmpts == 4
    
    MAP0(1,:)=[0 0 1];
    MAP0(2,:)=[0 1 0];
    MAP0(3,:)=[0 1 1];
    MAP0(4,:)=[1 0 0];
    
elseif nbcmpts == 5
    
    MAP0(1,:)=[0 0 1];
    MAP0(2,:)=[0 1 1];
    MAP0(3,:)=[0 1 0];
    MAP0(4,:)=[1 0 1];
    MAP0(5,:)=[1 0 0];
    
elseif nbcmpts == 6
    
    MAP0(1,:)=[0 0 1];
    MAP0(2,:)=[0 1 1];
    MAP0(3,:)=[0 1 0];
    MAP0(4,:)=[1 0 1];
    MAP0(5,:)=[1 1 0];
    MAP0(6,:)=[1 0 0];
   
else
    MAP0(1,:)=[0 0 1];
    MAP0(2,:)=[0 1 1];
    MAP0(3,:)=[0 1 0];
    MAP0(4,:)=[1 0 1];
    MAP0(5,:)=[1 1 0];
    MAP0(6,:)=[1 0 0];
    for j=7:nbcmpts
        tempmapj=rand(1,6)*MAP0(1:6,:);
        MAP0(j,:)=tempmapj-floor(tempmapj);
    end
end


% MAP0=MAP1-MAP0;
