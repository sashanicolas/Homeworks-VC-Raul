clc;
close all;

%carrega imagens
img1 = imread('img/rsz_imga.png');
img2 = imread('img/rsz_imgb.png');

%roda o SURF
% Get the Key Points
  Options.upright=true;
  Options.tresh=0.0001;
  Ipts1=OpenSurf(img1,Options);  
  Ipts2=OpenSurf(img2,Options);
% Put the landmark descriptors in a matrix
  D1 = reshape([Ipts1.descriptor],64,[]); 
  D2 = reshape([Ipts2.descriptor],64,[]); 
% Find the best matches
  err=zeros(1,length(Ipts1));
  cor1=1:length(Ipts1); 
  cor2=zeros(1,length(Ipts1));
  
  for i=1:length(Ipts1),
      distance=sum((D2-repmat(D1(:,i),[1 length(Ipts2)])).^2, 1);
      [err(i),cor2(i)]=min(distance);
          %  min(distance)
  end
% Sort matches on vector distance
  [err, ind]=sort(err); 
  cor1=cor1(ind); 
  cor2=cor2(ind);
% Show both images
  I = zeros([size(img1,1) size(img1,2)*2 size(img1,3)]);
  I(:,1:size(img1,2),:)=img1; I(:,size(img1,2)+1:size(img1,2)+size(img2,2),:)=img2;
  figure, imshow(I/255); hold on;
% Show the best matches
  for i=1:30,
      c=rand(1,3);
      plot([Ipts1(cor1(i)).x Ipts2(cor2(i)).x+size(img1,2)],[Ipts1(cor1(i)).y Ipts2(cor2(i)).y],'-','Color',c);
      plot([Ipts1(cor1(i)).x Ipts2(cor2(i)).x+size(img1,2)],[Ipts1(cor1(i)).y Ipts2(cor2(i)).y],'o','Color',c);
  end
%trata e pega pegar os pontos do SURF, bota em
x2Trans = [];
y2Trans = [];
xBase = [];
yBase = [];
for i=1:30,
    x2Trans = [x2Trans; Ipts1(cor1(i)).x];
    y2Trans = [y2Trans; Ipts1(cor1(i)).y];
end
for i=1:30,
    xBase = [xBase; Ipts2(cor2(i)).x];
    yBase = [yBase; Ipts2(cor2(i)).y];
end

%x2Trans
%y2Trans
%xBase
%yBase

%calcula H com uma funcao ransac
hTeste = ransac(x2Trans,y2Trans,xBase,yBase, 20, 100, 1, 10, length(x2Trans)); %6, 1000, 0.3, 4
                                             %n  %k  %t  %d  %s
hTeste
%gera a panoramica com a funcao Panorama2
[imgPanteste, imgOut1, imgOut2] = Panorama2(img1,img2,hTeste);

imwrite(imgPanteste,'img/panSURF.png');
%%}

