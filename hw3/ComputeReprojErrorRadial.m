function [Erro,pReproj]=ComputeReprojErrorRadial(P,p,M,k1,k2,k3)
pReproj=p;
for i=1:size(P,2)
    % Computes the squared distance from the image center
    a1=M(1,1:3)';
    a2=M(2,1:3)';

    d2=p(1,i)*a2-p(2,i)*a1;
    d2=d2'*d2;
    a1xa2=cross(a1,a2);
    d2=d2/(a1xa2'*a1xa2);
    
    % Computes the reprojection error
    lambda=1+k1*d2+k2*d2^2+k3*d2^3;
    pReproj(:,i)=diag([1/lambda 1/lambda 1])*M*P(:,i);
    pReproj(:,i)=pReproj(:,i)/pReproj(3,i);
end
% [k1 k2 k3 ]
Erro=mean(sqrt(sum((pReproj-p).^2)));

return

    