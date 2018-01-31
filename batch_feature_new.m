function [ succ]=batch_feature_new(input_folder,output_filename,dt,Tmax)
% This code is the implementation of the 3D steganalytic feature set,
% LFS76, proposed in "Steganalysis of 3D Objects Using Statistics of Local Feature Sets"
% Information Sciencecs, Volumes 415¨C416, November 2017, Pages 85-99
% by the authors of Zhenyu Li and Adrian G. Bors
%
% If you have any questions about the code, please contact
% zheenyuli@gmail.com

%INPUT:
%input_folder - the path of all the 3D objects in .off format
%output_filename - the path of the extracted feature
%dt - parameter for the speed of Laplacian smoothing, corresponding to
%\lambda in the paper.
%Tmax - parameter for the time of Laplacian smoothing, number of
%iteration=Tmax/dt.
%OUTPUT:
%succ - number features successfully extracted from the object 

[tmp, lghPath] = size(input_folder);
if input_folder(lghPath) ~= '\'
    input_folder = strcat(input_folder, '\');
end



drt = dir(input_folder);
[fnum, tmp] = size(drt);

FEA=zeros(fnum-2,76);
Fname=cell(fnum-2,1);
succ = 0;

for i = 1:fnum
    if drt(i).isdir == 1
        continue;
    end
    [tmp,lghName] = size(drt(i).name);
     if (isempty(findstr(lower(drt(i).name), '.off'))) || ((findstr(lower(drt(i).name), '.off') + 3) ~= lghName)
            continue;
     end
    inName = strcat(input_folder, drt(i).name);
    ori_file=inName;
    [FEA(i-2,:)] = LFS76_fea( ori_file,dt,Tmax ); 
 
    Fname{i-2} = drt(i).name;
    save(output_filename,'FEA', 'Fname','-v7.3');
    succ=succ+1
 
end
end
