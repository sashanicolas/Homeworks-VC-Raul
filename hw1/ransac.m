function H = ransac(u2Trans,v2Trans,uBase,vBase)
% H=ransac(u2Trans,v2Trans,uBase,vBase,)
% Computes the homography H applying the Direct Linear Transformation
%
% INPUTS: 
% u2Trans, v2Trans - vectors with coordinates u and v of the image to be transformed (p')
% uBase, vBase - vectors with coordinates u and v of the base image p 
%
% OUTPUT
% H - 3x3 matrix with the Homography
%
% Sasha Nicolas - 16/09/2014
    
    n = 6; %menor numero de pontos necessarios
    k = 3; %numero de iteracoes necessarias
    t = 100; %limite para decidir se ponto fits well
    d = 4; %numero de pontos proximos necessarios
    s = 18; %numero de pontos totais
    maiorQntInliers = 0;
    H=[];
    
    uns = repmat(1, 1, n);
    uns2 = repmat(1, 1, s-n);
    for i = 1:1
        %pega n random pontos
        indB = randIndex(s,n);
        indB = [1 2 3 4 5 6];
        indT = indB; %randIndex(s,n);
        %indT = [1 2 3 4 5 6];
        %calculca o h para os pontos randomicos
        h = DLTnorm(uBase(indB),vBase(indB),u2Trans(indT),v2Trans(indT));
        %calcula a distancia para cada ponto usando o h
        pontosB = [uBase(indB)'; vBase(indB)'];
        %pontosB
        %h
        pontosT = h*[u2Trans(indT)'; v2Trans(indT)';uns];
        %pontosT
        %divide pelo w da coord homogenea
        pontosT = bsxfun (@rdivide, pontosT, pontosT(3,:)); 
        %pontosT
        distRandom = sqrt( (pontosB(1,1:end)-pontosT(1,1:end)).^2 + (pontosB(2,1:end)-pontosT(2,1:end)).^2 );
        distRandom
        
        %para cada ponto fora da amostragem, calcula a distancia com o h
        %determinado
        indFora = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
        indFora(indB) = [];
        indFora
        pontosForaB = [uBase(indFora)'; vBase(indFora)'];
        pontosForaT = h*[u2Trans(indFora)'; v2Trans(indFora)';uns2];
        pontosForaT = bsxfun (@rdivide, pontosForaT, pontosForaT(3,:)); 
        distFora = sqrt( (pontosForaB(1,1:end)-pontosForaT(1,1:end)).^2 + (pontosForaB(2,1:end)-pontosForaT(2,1:end)).^2 );
        distFora
        indMaiores = find(distFora<t); %quant de dist menores que t
        qntMaiores = length(indMaiores);
        if qntMaiores >= d
            if qntMaiores > maiorQntInliers 
                %salva o H
                H = h;
                %junta os inliers indMaiores com os randomicos e calcula H
                %de novo
                
            end
        end
    end
    

return

function [uSelec, vSelec, uNsel, vNsel] = pegaRandomPontos(n, uPontos,vPontos)
    
 return
