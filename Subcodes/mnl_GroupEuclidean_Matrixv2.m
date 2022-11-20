function [EuD_Matrix]=mnl_GroupEuclidean_Matrixv2(data1,data2)
%Atempt to speed up the calculation of Euclidean distances
%% New Method
sz1=size(data1,1);
sz2=size(data2,1);
Dist_Matrix=nan(sz1,sz2);
%Calculate the distances
for i=1:sz1
    temp=data2-data1(i,:); %Sum
    temp=temp.^2; %Square it to remove negatives
    a=sum(temp,2); %sum of squares
    Dist_Matrix(i,:)=a'; %Invert and add to distance matrix
    clear temp a
end
EuD_Matrix=Dist_Matrix.^0.5;
end