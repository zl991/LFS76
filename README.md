This code is the implementation of the 3D steganalytic feature set,
LFS76, proposed in "Steganalysis of 3D Objects Using Statistics of Local Feature Sets"
Information Sciencecs, Volumes 415–416, November 2017, Pages 85-99
by the authors of Zhenyu Li and Adrian G. Bors

The code was implemented by Zhenyu Li, a PhD student from University of York, 
under the supervison of Dr. Adrian G. Bors (https://github.com/AdrianBors).
If you have any questions about the code, please contact
zheenyuli@gmail.com

batch_feature_new.m is the batch running file.
LFS76_fea.m is the main function extracting the 76 dim local feature set.
Many functions that used in the codes from other researchers' work and they are listed
in the following.

The functions that calculate the face normals, vertex normals and curvatures are implemented by Itzik Ben Shabat, see 
https://cn.mathworks.com/matlabcentral/fileexchange/47134-curvature-estimationl-on-triangle-mesh?s_tid=prof_contriblnk

Some functions are from the Toolbox Graph implemented by Gabriel Peyre, see
https://cn.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph?s_tid=prof_contriblnk

Some functions are from the geom3d implemented by  David Legland, see
https://cn.mathworks.com/matlabcentral/fileexchange/24484-geom3d?s_tid=prof_contriblnk

Author Zhenyu Li
