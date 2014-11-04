function [alfa,beta,theta,u0,v0,M,K,R,t,Kappa]=CameraCalibrationDemo
% [alfa,beta,theta,u0,v0,M,K,R,t,Kappa]=CameraCalibrationDemo
% This is an iteractive tool that allows the user to enter control points in an
% image of a calibration rig, and returns the camera parameters. 
%
% OUTPUTS:
% alfa   - focal length in the horizontal direction in pixesls
% beta   - focal length in the vertical direction in pixesls
% theta  - skew in degrees
% u0     - the u coordinate of the origin in image frame 
% v0     - the v coordinate of the origin in image frame
% M      - Projection Matrix
% K      - Matrix of intrinsic Parameters
% R      - the rotation matrix
% t      - the translation vector
% Kappa  - 1x3 vector containing the coefficients of radial distortion
%
% Raul Queiroz Feitosa Junho/2009

%%

Radial=input('Consider radial distortion (y/n)? - ','s');

RigImage=input('Enter the name of the image file - ','s');
figure(1);
subplot(1,1,1);
imshow(RigImage);
 
XMAX=input('Enter the number of esquares in each plane of the Calibration Rig - ');
XMAX=XMAX-1;

%% Collects Points

P=[ 0 0 0;...
    0 0 1;...
    0 0 2;...
    0 1 0;...
    0 1 1;...
    0 1 2;...
    0 2 0;...
    0 2 1;...
    0 2 2];
P
P=[P ;P(:,[2 1 3]);P(:,[3 2 1])];
P
P=unique(P,'rows')/2;
P

Coleta=input('Do you want to use pre stored points (y/n)? - ','s');

captura=0;
if Coleta(1)=='y',
    NomeArquivo=input('.mat file containing the points (without extension) - ','s');
    try 
        eval(['load ' NomeArquivo]);
    catch
        captura=1;
    end
    if captura==1 || size(p,1)~=3 || size(p,2)~=size(P,2) || size(P,1)~=4,
        disp('Input file has inconsistent data - please select a new set of points')
        captura=1;
    else
        captura=0;
    end
else
    captura=1;
end
figure(1);
if captura==1,
    p=[];
    for pt=1:size(P,1);
        title(['Zoom in the point [x y z] = ['...
            num2str(floor(XMAX*P(pt,:))) '] and press the spacebar' ]);
        zoom on
        pause
        [u,v,B]=ginput(1);
        p=[p;u v 1];
        hold on
        plot(u,v,'r+');
        zoom out
    end
    title('Corresponding Points have been collected!')
    L=input('Enter the square side in cm - ');
    P=[L*floor(XMAX*P');ones(1,size(P,1))];

    p=p';
    p(2,:)=-p(2,:);
    save LastCollectedPoints P p
    disp ('Your points were saved in file "LastCollectedPoints.m"');
end

%% Computes the Projection Matrix without Radial Distortion
if Radial(1)~='y' %sem distorcao radial
    
    Kappa=[0 0 0];
    %% Compute the Projection matrix and extracts parameters
    M=GoldStandardCameraCalibration(P,p);
    [K,R,t]=krt(M);
    M=K*[R t];

    alfa=K(1,1);
    theta=acotd(-K(1,2)/alfa);
    if theta<0,
        theta=theta+180;
    end
    beta=K(2,2)*sind(theta);
    u0=K(1,3);
    v0=K(2,3);

%%  Reprojects the points
    pReproj=M*P;
    pReproj=pReproj./repmat(pReproj(3,:),[3 1]);
    Erro=(pReproj-p);
    Erro=Erro.^2;
    Erro=sum(Erro);
    Erro=sqrt(Erro);
    Erro=mean(Erro);
    
%% Computes the Projection Matrix with Radial Distortion
else
    [alfa,beta,theta,u0,v0,M,K,R,t,k1,k2,k3,pReproj,Erro]=NonLinearCalibration2(P,p);
    Kappa=[k1 k2 k3];
end
%%  Paints the Reprojected points
figure;
hold off
imshow(RigImage);
hold on
plot(p(1,:),-p(2,:),'r+');
hold on
plot(pReproj(1,:),-pReproj(2,:),'g+');
plot(u0,-v0,'yo')
legend('Captured','Reprojected','Image Center')
title(['Mean Reprojection Error = ' num2str( Erro,2) ' pixel - Radial Distortion = ' num2str(Kappa)])

return

    