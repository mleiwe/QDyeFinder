function [M]=mnl_PermsNoColumnBias(v,n)
% Function that works out the different combinations of selecting n values
% from vector v, regardless of order
%Inputs
% v - vector e.g. [1 2 3 4 5 6 7]
% n - number to select
%
%Outputs
% M - Matrix with all the possibilities

ap=perms(v); %All permeabilities
uap=unique(ap(:,1:n),'rows'); %Remove the ones we aren't interested in dimensions and also identical
%% Now for the tricky part - Removing combinations regardless of column
cc=1;
M=[];
for i=1:length(uap) %For each row
    Val1=uap(i,:);
    IdenticalMatrix=zeros(length(uap),n);
    for j=1:n %For each value of this row
        source=Val1(j);
        %Find where this value appears
        [r,c]=find(uap==source);
        for k=1:size(r,1)
            IdenticalMatrix(r(k),c(k))=1;
        end
    end
    %Now identify the identical ones (the matrix should add up to n)
    IdenticalMatrix(i,:)=zeros(1,n);%First blank the self comparison - Not sure if necessary
    Test=sum(IdenticalMatrix,2); %Sum along the columns
    index=find(Test==n);
    %Is this the first one
    if i<=min(index)
        M(cc,:)=Val1;
        cc=cc+1;
    end
    clear index Test
end
end

