function H = ransac(u2Trans,v2Trans,uBase,vBase, this_n, this_k, this_t, this_d, this_s)
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
    
    n = this_n; %6; %menor numero de pontos necessarios
    k = this_k; %3; %numero de iteracoes necessarias
    t = this_t; %100; %limite para decidir se ponto fits well
    d = this_d; %4; %numero de pontos proximos necessarios
    s = this_s; %18; %numero de pontos totais
    maiorQntInliers = 0;
    H=[];
    
    %k=1;
    %n=11;
    
    uns = repmat(1, 1, n);
    uns2 = repmat(1, 1, s-n);
    menor = 1000000;
    for i = 1:k
        %pega n random pontos
        indB = randIndex(s,n);
        %indB = [5, 11, 13, 17, 18, 19, 20, 24, 33, 34, 36];
        %calculca o h para os pontos randomicos
        h = DLTnorm(uBase(indB),vBase(indB),u2Trans(indB),v2Trans(indB));
        %calcula a distancia para cada ponto usando o h
        pontosB = [uBase(indB)'; vBase(indB)'];
        pontosT = h*[u2Trans(indB)'; v2Trans(indB)';uns];
        %divide pelo w da coord homogenea
        pontosT = bsxfun (@rdivide, pontosT, pontosT(3,:)); 
        distRandom = sqrt( (pontosB(1,1:end)-pontosT(1,1:end)).^2 + (pontosB(2,1:end)-pontosT(2,1:end)).^2 );
        %distRandom
        if min(distRandom)<menor
           menor = min(distRandom);
        end
        
        %para cada ponto fora da amostragem, calcula a distancia com o h
        %determinado
        indFora = [];
        for iterador = 1 : s
            indFora(iterador) = iterador;
        end
        indFora(indB) = [];
        %indFora
        pontosForaB = [uBase(indFora)'; vBase(indFora)'];
        pontosForaT = h*[u2Trans(indFora)'; v2Trans(indFora)';uns2];
        pontosForaT = bsxfun (@rdivide, pontosForaT, pontosForaT(3,:)); 
        distFora = sqrt( (pontosForaB(1,1:end)-pontosForaT(1,1:end)).^2 + (pontosForaB(2,1:end)-pontosForaT(2,1:end)).^2 );
        %distFora
        indMaiores = find(distRandom<t); %quant de dist menores que t
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
    
    menor

return

