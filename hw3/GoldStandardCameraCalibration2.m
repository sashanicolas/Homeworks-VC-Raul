function M=GoldStandardCameraCalibration2(P,p)
% M=GoldStandardCameraCalibration(P,p)
% Computes the perspective projection matrix applying the Gold Standard
% algorithm.
%
% INPUTS:
% P - a 4xn matrix containing the calibration points in the world frame
% p - a 3xn matrix containing the corresponding points in the image frame
%
% OUTPUT
% M - the projection matrix
%
% Raul Queiroz Feitosa -Maio/2005
% Edit by Sasha Nicolas November/2014

% -------------------------------------------------------------------------
% Normalizing the image points
% -------------------------------------------------------------------------
[p,T]=NormalizePoints(p);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Normalizing the space points
% -------------------------------------------------------------------------
[P,U]=NormalizePoints(P);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Builds the matrix PP - que eh a matrix com os pontos e os zeros
% -------------------------------------------------------------------------
PP=[];
for ponto=1:size(P,2);
    PP=[PP; P(:,ponto)' 0 0 0 0 ];
    PP=[PP; 0 0 0 0 P(:,ponto)' ];
end

%arruma pzinho em b
b=[];
for ind=1:size(p,2);
    b=[b;p(1,ind);p(2,ind)];
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Usa  PseudoInversa
% -------------------------------------------------------------------------
Px = (inv(PP'*PP)*PP');
m=Px*b;

Mpseudo=[m([ 1 2 3 4;5 6 7 8]); 0 0 0 1];

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Applies Non Linear Refinement
% -------------------------------------------------------------------------
try
    options = optimset('MaxIter',1,...
    'Display','off',...
    'TolFun',1e-300,...
    'Algorithm', 'levenberg-marquardt',...
    'DiffMinChange',1e-40);
catch % for earlier MATLAB versions
    options = optimset('MaxIter',1,...
    'Display','off',...
    'TolFun',1e-300,...
    'DiffMinChange',1e-40);
end

[MnonLin,resnorm,residual,exitflag] = ...
    lsqnonlin(@(M) ReprojectionError(M,P,p),Mpseudo,[],[],options);
% -------------------------------------------------------------------------

if ~any(MnonLin(:)~=Mpseudo(:))
    disp('The nonlinear method brought no improvement in relation to the DLT!');
end

% -------------------------------------------------------------------------
% Denormalization
% -------------------------------------------------------------------------
M=inv(T)*MnonLin*U;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Makes M(3,1:3) a unitary vector
% -------------------------------------------------------------------------
M=M/sqrt(sum(M(3,1:3).^2));
% -------------------------------------------------------------------------

return


 function F = ReprojectionError(M,P,p)
        
     pfromM=M*P;
     pfromM=pfromM./repmat(pfromM(3,:),[3,1]);
     F=p-pfromM;
     F=sqrt(sum(F(1:2,:).^2));

return



