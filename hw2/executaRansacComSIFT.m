clc;
close all;

%imagem base?
%img1 = imread('img/imga.jpg');
%img1 = imread('img/rsz_imga.png');
img1 = imread('img/rsz_pan_a_b_2.png');

%imagem a ser trasnformada?
%img2 = imread('img/imgb.jpg');
img2 = imread('img/rsz_imgc.png');

%roda o sift
%[num,loc1,loc2] = matchSIFTRaul('img/rsz_pan_a_b_2.png', 'img/rsz_imgc.png');

%trata e pega pegar os pontos do sift, bota em
y2Trans = loc1(:,1);
y2Trans
x2Trans = loc1(:,2);
yBase = loc2(:,1);
xBase = loc2(:,2);

%calcula H com uma funcao ransac
hTeste = ransac(x2Trans,y2Trans,xBase,yBase, 14, 100, 1, 10, length(loc1)); %6, 1000, 0.3, 4
                                             %n  %k   %t   %d  %s
hTeste
%gera a panoramica com a funcao Panorama2
[imgPanteste, imgOut1, imgOut2] = Panorama2(img1,img2,hTeste);

%imwrite(imgPanteste,'img/panSIFT.png');
%%}
