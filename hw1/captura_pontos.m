function [x2Trans y2Trans xBase yBase] = captura_pontos(I2Trans,IBase)
% function [x2Trans y2Trans xBase yBase] = captura_pontos(I2Trans,IBase)
% Captura pontos de 2 imagens estereo manualmente
% 
% Entradas da Funcao:
% I2Trans - imagem da esquerda
% IBase - imagem da direita
%
% Saídas da Funcao:
% x2Trans - vetor com as coordenadas x dos pontos, na imagem da esquerda
% y2Trans - vetor com as coordenadas y dos pontos, na imagem da esquerda
% xBase - vetor com as coordenadas x dos pontos, na imagem da direita
% yBase - vetor com as coordenadas y dos pontos, na imagem da direita
%
% O arquivo "pontos_homologos" conterá as saídas da função
%
% obs: Para terminar a captura dos pontos, clicar com o botão direito!
%
% Raul Queiroz Feitosa - Março/2010

x2Trans= [];
y2Trans = [];
xBase = [];
yBase = [];
figure(1)
h1 = subplot(1,2,1);
imshow(IBase)
h2 = subplot(1,2,2);
imshow(I2Trans)

but=1;

while but==1 || length(xBase)<4    
    axes(h1)
    title('Faça zoom na área de interesse e tecle "ENTER"');
    pause
    
    axes(h1)
    title({'Imagem BASE';'Selecione um ponto'})
    hold on
    [xi, yi, but] = ginput(1);
    xBase = [xBase;xi];
    yBase = [yBase;yi];
    hold on
    plot(xi,yi,'r+');
    title('')
    
    axes(h2)  
    title({'Imagem 2Trans - Selecione o homólogo';'para terminar botão direito'})
    hold on
    [xi, yi, but] = ginput(1);
    x2Trans = [x2Trans;xi];
    y2Trans = [y2Trans;yi];
    hold on
    plot(xi,yi,'g+');
    title('')
end

save pontos_homologos2 x2Trans y2Trans xBase yBase

return
