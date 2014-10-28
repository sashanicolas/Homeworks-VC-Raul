function H = DLTnorm(u2Trans,v2Trans,uBase,vBase)
% H=DLT(u2Trans,v2Trans,uBase,vBase,)
% Computes the homography H applying the Direct Linear Transformation
% The transformation is such that
% p = H p'
% (uBase vBase 1)'=H*(u2Trans v2Trans 1)'
%
% INPUTS: 
% u2Trans, v2Trans - vectors with coordinates u and v of the image to be transformed (p')
% uBase, vBase - vectors with coordinates u and v of the base image p 
%
% OUTPUT
% H - 3x3 matrix with the Homography
%
% Sasha Nicolas - 16/09/2014
    
    %normaliza
    [uBase, vBase, T1] = normaliza(uBase, vBase);
    [u2Trans, v2Trans, T2] = normaliza(u2Trans, v2Trans);
    
    qntPontos = length(u2Trans); %numero de pontos, tanto faz pegar de qq um
    
    M = zeros(3*qntPontos,9); %inicializa matriz com zeros
    
    zero = [0 0 0]; %vetor de zeros
    for n = 1:qntPontos
        %pega o ponto n em imgBase
        %X = x1(:,n)';
        X = [uBase(n) vBase(n) 1];
        %X
        %pega o ponto n em imgTrans
        x = u2Trans(n); %x2(1,n); 
        y = v2Trans(n); %x2(2,n); 
        w = 1; %x2(3,n);
        
        M(3*n-2,:) = [  zero  -w*X  y*X];
        M(3*n-1,:) = [ w*X   zero  -x*X];
        M(3*n  ,:) = [-y*X  x*X   zero ];
    end
    
    %resolver o sistema
    [U,D,V] = svd(M,0); % 'Economy' decomposition for speed
    % Extract homography
    H = reshape(V(:,9),3,3)';
    % Denormalise
    H = T2\H*T1;
    
return

function [uNovos, vNovos, T] = normaliza(uPontos,vPontos)
    uPontos=uPontos';
    vPontos=vPontos';
    %media de x
    m1 = mean(uPontos(1,1:end)); 
   
    %media de y
    m2 = mean(vPontos(1,1:end)); 
   
    %move todos
    uAux = uPontos(1,1:end) - m1;
    vAux = vPontos(1,1:end) - m2;
    %distancia
    dist = sqrt(uAux.^2 + vAux.^2);
    meandist = mean(dist(:));  % Ensure dist is a column vector for Octave 3.0.1
    
    scale = sqrt(2)/meandist;
    
    T = [scale   0   -scale*m1
         0     scale -scale*m2
         0       0      1      ];
    
    uns = ones(1, length(uPontos));
    
    mat = [uPontos(1,1:end); vPontos(1,1:end); uns];
    
    mat = T*mat;
    %mat
    uNovos = mat(1,:);
    vNovos = mat(2,:);
return