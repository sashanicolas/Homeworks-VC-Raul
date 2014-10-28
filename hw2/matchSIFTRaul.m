% [num,loc1,loc2] = matchSIFTRaul(image1, image2)
%
% This function reads two images, finds their SIFT features, and
%   displays lines connecting the matched keypoints.  A match is accepted
%   only if its distance is less than distRatio times the distance to the
%   second closest match.
% It returns the number of matches displayed, and the coordinates of the
% matched points in both images (loc1,loc2)
%
% Example: match('scene.pgm','book.pgm');

function [num,loc1,loc2] = matchSIFTRAUL(image1, image2)

% Find SIFT keypoints for each image
[im1, des1, loc1] = sift(image1);
[im2, des2, loc2] = sift(image2);

% For efficiency in Matlab, it is cheaper to compute dot products between
%  unit vectors rather than Euclidean distances.  Note that the ratio of 
%  angles (acos of dot products of unit vectors) is a close approximation
%  to the ratio of Euclidean distances for small angles.
%
% distRatio: Only keep matches in which the ratio of vector angles from the
%   nearest to second nearest neighbor is less than distRatio.
distRatio = 0.6;   

% For each descriptor in the first image, select its match to second image.
des2t = des2';                          % Precompute matrix transpose
for i = 1 : size(des1,1)
   dotprods = des1(i,:) * des2t;        % Computes vector of dot products
   [vals,indx] = sort(acos(dotprods));  % Take inverse cosine and sort results

   % Check if nearest neighbor has angle less than distRatio times 2nd.
   if (vals(1) < distRatio * vals(2))
      match(i) = indx(1);
   else
      match(i) = 0;
   end
end

% Commented by Raul<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% Create a new image showing the two images side by side.
im3 = appendimages(im1,im2);

% Show a figure with lines joining the accepted matches.
figure('Position', [100 100 size(im3,2) size(im3,1)]);
colormap('gray');
imagesc(im3);
hold on;
cols1 = size(im1,2);

Hpar=gcf;
% figure;imshow(im1);
% H1=gcf;
% figure;imshow(im2);
% H2=gcf;

x1=[];
y1=[];
x2=[];
y2=[];
for i = 1: size(des1,1)
  if (match(i) > 0)
      figure(Hpar)
      line([loc1(i,2) loc2(match(i),2)+cols1], ...
          [loc1(i,1) loc2(match(i),1)], 'Color', 'c');
      
      x1=[x1 loc1(i,2)];
      y1=[y1 loc1(i,1)];  
      x2=[x2 loc2(match(i),2)];
      y2=[y2 loc2(match(i),1)];
      
%       figure(H1);
%       hold on    
%       plot(x1(end),y1(end),'r+');      
%       figure(H2)
%       hold on
%       plot(x2(end),y2(end),'g+');
          
  end
end
axis off;
hold off;

indices=find(match>0);
loc1=[y1(:) x1(:)];
loc2=[y2(:) x2(:)];
num = length(loc1);
fprintf('Found %d matches.\n', num);




