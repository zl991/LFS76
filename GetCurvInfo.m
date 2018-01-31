function [ GausianCurvature, CurvatureRatio] = GetCurvInfo( vertex,face, getderivatives )
%GETCURVINFO Summary of this function goes here
%   Detailed explanation goes here
    FV.faces = face';
    FV.vertices=vertex';
    
    [PrincipalCurvatures,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude]= GetCurvatures( FV , getderivatives);
    
    GausianCurvature=PrincipalCurvatures(1,:).*PrincipalCurvatures(2,:);
    
    CurvatureRatio=min(abs(PrincipalCurvatures))./max(abs(PrincipalCurvatures));

end

