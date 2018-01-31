function [ F ] = LFS76_fea( ori_file,dt,Tmax )
% This code is the implementation of the 3D steganalytic feature set,
% LFS76, proposed in "Steganalysis of 3D Objects Using Statistics of Local Feature Sets"
% Information Sciencecs, Volumes 415¨C416, November 2017, Pages 85-99
% by the authors of Zhenyu Li and Adrian G. Bors
%
% If you have any questions about the code, please contact
% zheenyuli@gmail.com

%INPUT:
%ori_file - the path of the 3D object in .off format
%dt - parameter for the speed of Laplacian smoothing, corresponding to
%\lambda in the paper.
%Tmax - parameter for the time of Laplacian smoothing, number of
%iteration=Tmax/dt.
%OUTPUT:
%F - a 1X76 matrix of the local feature set

% read in the 3D object 
[vertex_o, face_o] = read_off(ori_file);
 
%%Laplacian smoothing the object 
options.symmetrize = 0;
options.normalize = 1;

% laplacian_type = 'distance';
laplacian_type = 'combinatorial';
% laplacian_type = 'conformal';

%  L can be a laplacian or a string containing the type of laplacian.
disp('--> Computing laplacian');
L = compute_mesh_laplacian(vertex_o,face_o,laplacian_type,options);

%dt represents the \lambda in the paper, which controls the speed of the
%diffusion
options.dt = dt;
%Tmax represents the time of the diffusion
options.Tmax = Tmax; 

%vertex_s represent the coordinates of the vertex in the smoothed object
vertex_s = perform_mesh_heat_diffusion(vertex_o,face_o,L,options);
face_s=face_o;

%Do the roation and normalization to the original object and its smoothed
%version
vertex_o=preprocess(vertex_o);
vertex_s=preprocess(vertex_s);


 

%calculate the adjacence graph
W = triangulation2adjacency(face_o);
 
 
%vertex coordinates in Cartesian coordiante system 
diff_vertex_position = abs(vertex_o-vertex_s)+eps;

%vertex norm in Cartesian coordiante system 
diff_vertex_longth= abs(sqrt(sum(vertex_o(:,:).^2))-sqrt(sum(vertex_s(:,:).^2)))+eps;

%Transform the object from Cartesian coordinate system to the Laplacian 
%coordinate system
L_o = (diag(sum(W,2)) - W)*vertex_o';
L_o=L_o';
L_s = (diag(sum(W,2)) - W)*vertex_s';
L_s=L_s';
%vertex coordinates in Laplacian coordiante system 
Ldiff_vertex_position = abs(L_o(:,:)-L_s(:,:))+eps;
%vertex coordinates in Laplacian coordiante system 
Ldiff_vertex_longth = abs(sqrt(sum(L_o(:,:).^2))-sqrt(sum(L_s(:,:).^2)))+eps;

%Angles between face normals of the object and its smoothed version
normals_o = faceNormal(vertex_o', face_o');
normals_s = faceNormal(vertex_s', face_s');
diff_faces = vectorAngle3d(normals_o,normals_s)+eps;
diff_faces(isnan(diff_faces))=eps;


FV_o.faces=face_o';
FV_s.faces=face_s';
FV_o.vertices=vertex_o';
FV_s.vertices=vertex_s';
%Calculate the vertex normals
[VertexNormals_o,~,~,~,~]=CalcVertexNormals(FV_o,normals_o);
[VertexNormals_s,~,~,~,~]=CalcVertexNormals(FV_s,normals_s);

%Angles between the vertex normals 
diff_vertex_normals = vectorAngle3d(VertexNormals_o,VertexNormals_s)+eps;
diff_vertex_normals(isnan(diff_vertex_normals))=eps;

%Get curvature information
getderivatives=0;
[GaussianCurvature_o, CurvatureRatio_o]=GetCurvInfo(vertex_o,face_o,getderivatives);
[GaussianCurvature_s, CurvatureRatio_s]=GetCurvInfo(vertex_s,face_s,getderivatives);
%Calculate the Gaussian curvature
diff_GaussianCurv=abs(GaussianCurvature_o-GaussianCurvature_s)+eps;
%Calculate the curvature ratio
diff_CurvatureRatio=abs(CurvatureRatio_o-CurvatureRatio_s)+eps;
diff_GaussianCurv(isnan(diff_GaussianCurv))=eps;
diff_CurvatureRatio(isnan(diff_CurvatureRatio))=eps;

%Get the dihedral angles  
edge_o=meshEdges(face_o');
edge_s=edge_o; 
[dihedral_o] = meshDihedralAngles(vertex_o', edge_o, face_o');
[dihedral_s] = meshDihedralAngles(vertex_s', edge_s, face_s');
diff_dihedral = abs(dihedral_o-dihedral_s)+eps;
diff_dihedral(isnan(diff_dihedral))=eps;
 

%spherical coordinate features
origin_o = mean(vertex_o,2);
nvertex=size(vertex_o,2);
origin_o=origin_o*ones(1,nvertex);
%Transform the vertex into spherical coordinate system
[az_o,elev_o,rho_o] = enu2aer(vertex_o(1,:)-origin_o(1,:),vertex_o(2,:)-origin_o(2,:),vertex_o(3,:)-origin_o(3,:));
%Get the edge information
edge_o=meshEdges(face_o');
%Get the length of the edge
EAaz_o=abs(az_o(edge_o(:,1))-az_o(edge_o(:,2)));
EAel_o=abs(elev_o(edge_o(:,1))-elev_o(edge_o(:,2)));
Erho_o=abs(rho_o(edge_o(:,1))-rho_o(edge_o(:,2)));

origin_s = mean(vertex_s,2);
nvertex=size(vertex_s,2);
origin_s=origin_s*ones(1,nvertex);
%Transform the vertex into spherical coordinate system
[az_s,elev_s,rho_s] = enu2aer(vertex_s(1,:)-origin_s(1,:),vertex_s(2,:)-origin_s(2,:),vertex_s(3,:)-origin_s(3,:));
%Get the edge information
edge_s=meshEdges(face_s');
%Get the length of the edge
EAaz_s=abs(az_s(edge_s(:,1))-az_s(edge_s(:,2)));
EAel_s=abs(elev_s(edge_s(:,1))-elev_s(edge_s(:,2)));
Erho_s=abs(rho_s(edge_s(:,1))-rho_s(edge_s(:,2)));

%The differences between the vertex of the original mesh and its smoothed version in the SCS
diff_az=abs(az_o-az_s)+eps;
diff_elev=abs(elev_o-elev_s)+eps;
diff_rho=abs(rho_o-rho_s)+eps;

%The differences between the edge length the original mesh and its smoothed version in the SCS
diff_EAaz=abs(EAaz_o-EAaz_s)+eps;
diff_EAel=abs(EAel_o-EAel_s)+eps;
diff_Erho=abs(Erho_o-Erho_s)+eps;

%Calculate the first four statistical moments of the differences of the
%geometric features
F=zeros(1,76);

%vertex coordinates Cartesian
F(1,1:12) = [mean(log(diff_vertex_position')),var(log(diff_vertex_position')),...
    skewness(log(diff_vertex_position')), kurtosis(log(diff_vertex_position'))];

%vertex norm  Cartesian
F(1,13:16) = [mean(log(diff_vertex_longth')),var(log(diff_vertex_longth')),...
    skewness(log(diff_vertex_longth')), kurtosis(log(diff_vertex_longth'))];

%vertex coordinates Laplacian
F(1,17:28) = [mean(log(Ldiff_vertex_position')),var(log(Ldiff_vertex_position')),...
    skewness(log(Ldiff_vertex_position')), kurtosis(log(Ldiff_vertex_position'))];

%vertex norm Laplacian
F(1,29:32) = [mean(log(Ldiff_vertex_longth')),var(log(Ldiff_vertex_longth')),...
    skewness(log(Ldiff_vertex_longth')), kurtosis(log(Ldiff_vertex_longth'))];

%Face normal
F(1,33:36) = [mean(log(abs(diff_faces'))),var(log(abs(diff_faces'))),...
    skewness(log(abs(diff_faces'))), kurtosis(log(abs(diff_faces')))];
%Dihedral angles
F(1,37:40) = [mean(log(diff_dihedral')),var(log(diff_dihedral')),...
    skewness(log(diff_dihedral')), kurtosis(log(diff_dihedral'))];
%Vertex normals
F(1,41:44) = [mean(log(diff_vertex_normals')),var(log(diff_vertex_normals')),...
    skewness(log(diff_vertex_normals')), kurtosis(log(diff_vertex_normals'))];
%Gaussian curvauture
F(1,45:48) = [mean(log(diff_GaussianCurv')),var(log(diff_GaussianCurv')),...
    skewness(log(diff_GaussianCurv')), kurtosis(log(diff_GaussianCurv'))];
%Curvature ratio
F(1,49:52) = [mean(log(diff_CurvatureRatio')),var(log(diff_CurvatureRatio')),...
    skewness(log(diff_CurvatureRatio')), kurtosis(log(diff_CurvatureRatio'))];
%Vertex coordinate in the spherical coordinate system
F(1,53:56) = [mean(log(diff_az')),var(log(diff_az')),...
    skewness(log(diff_az')), kurtosis(log(diff_az'))];
F(1,57:60) = [mean(log(diff_elev')),var(log(diff_elev')),...
    skewness(log(diff_elev')), kurtosis(log(diff_elev'))];
F(1,61:64) = [mean(log(diff_rho')),var(log(diff_rho')),...
    skewness(log(diff_rho')), kurtosis(log(diff_rho'))];

%Edge length in the spherical coordinate system
F(1,65:68) = [mean(log(diff_EAaz')),var(log(diff_EAaz')),...
    skewness(log(diff_EAaz')), kurtosis(log(diff_EAaz'))];
F(1,69:72) = [mean(log(diff_EAel')),var(log(diff_EAel')),...
    skewness(log(diff_EAel')), kurtosis(log(diff_EAel'))];
F(1,73:76) = [mean(log(diff_Erho')),var(log(diff_Erho')),...
    skewness(log(diff_Erho')), kurtosis(log(diff_Erho'))];


end

