function H = DLT(u2Trans,v2Trans,uBase,vBase)
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
    
    qntPontos = length(u2Trans); %numero de pontos, tanto faz pegar de qq um
    %qntPontos 
    M = zeros(3*qntPontos,9); %inicializa matriz com zeros
    
    zero = [0 0 0]; %vetor de zeros
    for n = 1:qntPontos-1
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
    
return