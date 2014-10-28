function [Pnormalized,T]=NormalizePoints(Pin)
%[Pnormalized,T]=NormalizePoints(Pin)
% Transforms the points in Pin, by applying a transformation of the typ
%  PNormalized=T* Pin;
%
% The transformation T is computed so that the centroid of the normalized 
% points is at the origin and the average distance to the origin is equal 
% to the square root of the dimension of Pin
%
% INPUTS:
% Pin  - matrix containg on each column input points in homogeneous 
%        coordinates along the largest dimension
%
% OUTPUTS:
% Pnormlized - the normalized output Points in homogeneous coordinates
% T          - the transformation matrix
%
% Raul Queiroz Feitosa - April/2010



%=============== Check input data consistency  ==============

[M,N]=size(Pin);
if M>N ,
    error('The input points are not stored column-wise!');
end
if any(Pin(end,:)~=1),
    error('The input points are not in nonhomogeneous coordinates~');
end

%=============== Normalization   ==================
% Computes the normalization matrix
PinAvg=zeros(M-1,1);
for d=1:(M-1),
    PinAvg(d)=mean(Pin(d,:));
end

% Computes the average distance to centroid "PinAvgDist"
PinAvgDist=Pin-[repmat(PinAvg,[1,N]);ones(1,N)];
PinAvgDist=mean(sqrt(sum(PinAvgDist.^2)));

% Computes the scale factor "s"
s=sqrt(M-1)/PinAvgDist;

% Builds the transformation matrix
T=eye(M);
T(:,M)=[-PinAvg; 1/s];
T=s*T;

% Normalizes the Input points
Pnormalized=T*Pin;

end
    

% uav2Trans=mean(u);
% vav2Trans=mean(v);
% Dav2Trans=mean(((u-uav2Trans).^2+(v-vav2Trans).^2).^(1/2));
% s=sqrt(2)/Dav2Trans;
% T=s*[1 0 -uav2Trans;0 1 -vav2Trans;0 0 1/s];
% aux=[u(:)';v(:)'; ones(1,length(u))];
% aux= T*aux;
% uN=aux(1,:);
% vN=aux(2,:);
% return