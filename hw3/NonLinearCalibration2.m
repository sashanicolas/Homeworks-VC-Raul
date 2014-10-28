function [alfa,beta,theta,u0,v0,M,K,R,t,k1,k2,k3,pReproj,Erro]=NonLinearCalibration(P,p)
%[alfa,beta,theta,u0,v0,M,K,R,t,k1,k2,k3,pReproj,Erro]=NonLinearCalibration(P,p)
% extracts the camera parameters
% from a set of points (P) in world and their projection on the image (p)
%
% INPUTS:
% P - a 4xn matrix containing the calibration points in the world frame
% p - a 3xn matrix containing the corresponding points in the image frame
%
% OUTPUTS:
% alfa   - focal length in the horizontal direction in pixesls
% beta   - focal length in the vertical direction in pixesls
% theta  - skew in degrees
% u0     - the u coordinate of the origin in image frame 
% v0     - the v coordinate of the origin in image frame  
% ro     - scale factor  
% M      - Projection Matrix
% K      - Matrix of intrinsic Parameters
% R      - the rotation matrix
% t      - the translation vector
% Kappa  - vector containing the 3 coefficientes of the radial distortion
%
% Raul Queiroz Feitosa Abril/2006

%=============================================================
% Computes Linear Calibration to obtain initial solution
%=============================================================

% Computes the parameters using Linear Approach
Mlinear=GoldStandardCameraCalibration(P,p);
[Klin,Rlin,tlin]=krt(Mlinear);
M=Klin*[Rlin tlin];

u0=Klin(1,3);
v0=Klin(2,3);

p=p-repmat([u0 v0 0]',[1 size(p,2)]);

%=============================================================
% From this point on in this program image coordinates p have 
% the origin at the imagem center
%=============================================================

% Computes again the parameters with shifted origin
M=GoldStandardCameraCalibration(P,p);
[K,R,t]=krt(M);
M=K*[R t];

x0=[M(:)' 0 0 0];
xNorm=[M(:)' 1e-1 1e-2 1e-3];
%=============================================================
% From this point on in this program image coordinates p have 
% the origin at the imagem center
%=============================================================

% Call non linear equation system solver

% options = optimset('MaxIter', 500,...
%     'PlotFcns',@optimplotresnorm,...
%     'TolFun',1e-8,...
%     'MaxFunEvals',1000000 );

options = optimset('MaxIter', 500,...
    'TolFun',1e-8,...
    'MaxFunEvals',1000000);

[x, residuo] = lsqnonlin(@(x) ...
    SoErro(x,P,p,xNorm),...
    x0./xNorm,[],[],options);
[Erro,M,K,R,t,pReproj]=ErroNonLinear(x,P,p,xNorm);

% Prepara valores de saída
K(1,3)=u0;
K(2,3)=v0;
alfa=K(1,1);
theta=acotd(-K(1,2)/alfa);
beta=K(2,2)*sind(theta);
k1=x(end-2);
k2=x(end-1);
k3=x(end);

pReproj=pReproj+repmat([u0 v0 0]',[1 size(p,2)]);

return
%=============================================================

%=============================================================
% Objective function
%=============================================================
function Erro=SoErro(x,P,p,xNorm)
[Erro,M,K,R,t,pProj]=ErroNonLinear(x,P,p,xNorm);
return
%=============================================================

%=============================================================
% Function that computes the reprojection Error
%=============================================================
function [Erro,M,K,R,t,pProj]=ErroNonLinear(x,P,p,xNorm)
x=x.*xNorm;
M=zeros(3,4);
M(:)=x(1:12);
k1=x(13);
k2=x(14);
k3=x(15);

[K,R,t]=krt(M);
u0=K(1,3);
v0=K(2,3);
alfa=K(1,1);
theta=acotd(-K(1,2)/alfa);
beta=K(2,2)*sind(theta);
p=p-repmat([u0 v0 0]',[1 size(p,2)]);

u=p(1,:);
v=p(2,:);
d2=(u.^2)/(alfa^2) + (v.^2)/(beta^2) + (u.*v)*cosd(theta)/(alfa*beta);

lambda=1+k1*d2+k2*d2.^2+k3*d2.^3;
pProj=M*P;
for i=1:size(p,2);
    pProj(:,i)=diag([1/lambda(i) 1/lambda(i) 1])*pProj(:,i);
    pProj(:,i)=pProj(:,i)/pProj(3,i);
end

Erro=pProj-p;
Erro=sqrt(sum(Erro.^2));
Erro=mean(Erro);

return
%=============================================================