function [uN,vN,T]=NormalizaPontos(u,v)
%[uN,vN,T]=NormalizaPontos(u,v)
% Esta função aplica uma transformação nos pontos definidos em u e v, assim
%  [uN vN 1]'=T* [u v 1]'
%
% Inicialmente a transformação T é calculada de modo que a média dos pontos
% transformados esteja sobre a origem uN(medio)=vN(médio)=0, e a distância
% dos pontos transformados em relação à origem seja igual a sqrt(2)
%
% ENTRDAS:
% u - coordenadas x
% v - coordenadas y
%
% SAIDAS:
% uN - coordenadas x normalizadas
% vN - coordenadas y normalizadas
% T  - matriz contendo a transformação
%
% Raul Queiroz Feitosa - Março/2010


%=============== Normalização   ==================
% Calcula a matriz de normalização
uav2Trans=mean(u);
vav2Trans=mean(v);
Dav2Trans=mean(((u-uav2Trans).^2+(v-vav2Trans).^2).^(1/2));
s=sqrt(2)/Dav2Trans;
T=s*[1 0 -uav2Trans;0 1 -vav2Trans;0 0 1/s];
aux=[u(:)';v(:)'; ones(1,length(u))];
aux= T*aux;
uN=aux(1,:);
vN=aux(2,:);
return