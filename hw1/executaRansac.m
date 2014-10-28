clc;
%confusao de o que eh imagem base e qual vai ser transformada
%imagem base?
img1 = imread('img/rsz_imga.png');
%imagem a ser trasnformada?
img2 = imread('img/rsz_imgb.png');
%chama janela para capturar pontos
%captura_pontos(img2,img1);
%carrega o resultado em vetores
load('pontos_homologos_rsz_a_b.mat');

%%{
%adiciona outros pontos alem dos pontos coletados
x2Trans = [x2Trans; 900; 293; 147; 289; 720];
y2Trans = [y2Trans; 760; 627; 179; 817; 918];
xBase = [xBase; 92; 50; 147; 873; 53];
yBase = [yBase; 267; 282; 91; 93; 728];

%calcula H com uma funcao
hTeste = ransac(x2Trans,y2Trans,xBase,yBase);

%gera a panoramica com a funcao Panorama2
[imgPanteste, imgOut1, imgOut2] = Panorama2(img2,img1,hTeste);

imwrite(imgPanteste,'resultados/panRansac.png');

%%}
