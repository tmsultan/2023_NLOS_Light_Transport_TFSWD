% ----------------------------------------------------------------------- %
% makegrid() creates the 3D voxel space coordinates using leftbottomfront % 
% and righttopback (boundary points of the desired 3D cube). The grid size%
% is defined by voxelSize.                                                %
% The current version uses meshgrid (as opposed to nested loops) and      %
% returns the same output as the old version                              %
% ----------------------------------------------------------------------- %
function [voxelCoordinates, volumeSize] = makegrid(leftbottomfront, righttopback, voxelSize)

ax1 = leftbottomfront(1):voxelSize:righttopback(1);
ax2 = leftbottomfront(2):voxelSize:righttopback(2);
ax3 = leftbottomfront(3):voxelSize:righttopback(3);

[X, Y, Z] = meshgrid(ax2, ax3, ax1);

voxelCoordinates = [X(:), Y(:), Z(:)];
voxelCoordinates = [voxelCoordinates(:,3) voxelCoordinates(:,1) voxelCoordinates(:,2)];
volumeSize = [length(ax1) length(ax2) length(ax3)];

% The above code replaces the following code
% volumeSize=ceil(abs(rightbottomback-lefttopfront)./gridsize);
% vx=[1 0 0].*gridsize;
% vy=[0 1 0].*gridsize;
% vz=[0 0 1].*gridsize;
% npixels=volumeSize(1)*volumeSize(2)*volumeSize(3);
% voxelCoordinates=zeros(npixels,3);
% j=1;
% for x=1:volumeSize(1)
%     for y=1:volumeSize(2)
%         for z=1:volumeSize(3)
%             voxelCoordinates(j,:)=rightbottomback+(x-1).*vx+(y-1)*vy+(z-1)*vz;
%             j=j+1;
%         end
%     end
% end
% end
