function M=LinearCameraCalibration(P,p)
% M=LinearCameraCalibration(P,p)
% Computes the perspective projection matrix applying the linear method
%
% INPUTS:
% P - a 4xn matrix containing the calibration points in the world frame
% p - a 3xn matrix containing the corresponding points in the image frame
%
% OUTPUT
% M - the projection matrix
%
% Raul Queiroz Feitosa -Maio/2005


% Builds the matrix PP
PP=[];
for ponto=1:size(P,2);
    PP=[PP; P(:,ponto)' 0 0 0 0 -p(1,ponto)*P(:,ponto)'];
    PP=[PP; 0 0 0 0 P(:,ponto)' -p(2,ponto)*P(:,ponto)'];
end

[vet,val]=eig(PP'*PP);
[val,indices]=sort(diag(val));

m=vet(:,indices(1));

M=m([ 1 2 3 4;5 6 7 8; 9 10 11 12]);

return







