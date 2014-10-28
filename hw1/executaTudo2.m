clc;
%confusao de o que eh imagem base e qual vai ser transformada
%imagem base?
img1 = imread('img/IMG_20140911_090410.jpg');
%imagem a ser trasnformada?
img2 = imread('resultados/panRansac.png');
%chama janela para capturar pontos
%captura_pontos(img2,img1);
%carrega o resultado em vetores
load('pontos_homologos2.mat');

%%{
  
%prepara os valores para entrada da funcao homography2d
uns = repmat(1, 1, length(xBase));
x1 = [xBase';yBase';uns];
x2 = [x2Trans';y2Trans';uns];
x1 = x1(:,1:end-1);
x2 = x2(:,1:end-1);

%calcula H com uma funcao
%hTeste = homography2d(x1,x2);
%hTeste = DLT(x2Trans,y2Trans,xBase,yBase);
x2Trans = x2Trans(1:end-1);
y2Trans = y2Trans(1:end-1);
xBase = xBase(1:end-1);
yBase = yBase(1:end-1);
hTeste = DLTnorm(x2Trans,y2Trans,xBase,yBase);
hTeste
%gera a panoramica com a funcao Panorama2
[imgPanteste, imgOut1, imgOut2] = Panorama2(img1,img2,hTeste);
 
%%}
