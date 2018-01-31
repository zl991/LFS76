function [vertex_new] = preprocess( vertex )
% This code is the implementation of the 3D steganalytic feature set,
% LFS76, proposed in "Steganalysis of 3D Objects Using Statistics of Local Feature Sets"
% Information Sciencecs, Volumes 415¨C416, November 2017, Pages 85-99
% by the authors of Zhenyu Li and Adrian G. Bors
%
% If you have any questions about the code, please contact
% zheenyuli@gmail.com

%INPUT:
%vertex - a 3 X n matrix for the vertex coordinates of the 3D object, 
%n is the number of the vertex in the object  
%OUTPUT:
%vertex_new - a 3 X n matrix for the vertex coordinates of the 3D object 
%after the preprocessing

%transform the coordinate by PCA
trans=pca(vertex');
vertex_new=trans*vertex;

%scale the model in to unit cube 
[range_up]=max(vertex_new');
[range_low]=min(vertex_new');
range=range_up-range_low;
[range_max,axis]=max(range);

for i=1:3   
    vertex_new(i,:) = (vertex_new(i,:)-range_low(i)*ones(1,length(vertex_new)))/range_max;
end

end


