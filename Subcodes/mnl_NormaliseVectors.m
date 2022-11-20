function [Matrix]=mnl_NormaliseVectors(Matrix)
% Input
% Matrix - row of points,columns represent dimensions
% Output
% Matrix - same as above but now normalised
matrix_Norm=sqrt(sum((Matrix.*Matrix),2));
Matrix=Matrix./repmat(matrix_Norm,1,size(Matrix,2));
end
