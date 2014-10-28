function [Panorama, imBase, imTransfor]=Panorama2(ITransf,IBase,H)

T=maketform('projective',H');
[imTransf,xoffset,yoffset]=imtransform(ITransf,T);
[yB,xB]=size(IBase(:,:,1));
[yT,xT]=size(imTransf(:,:,1));
x=round(xoffset);
y=round(yoffset);

xini=min(x(1,1),1);
yini=min(y(1,1),1);
xfin=max(x(1,2),xB);
yfin=max(y(1,2),yB);

if xini>0
    length_panorama_x=xfin;
else
    length_panorama_x=xfin-xini+1;
end
if yini>0
    length_panorama_y=yfin;
else
    length_panorama_y=yfin-yini+1;
end

img_panorama=uint8(zeros(length_panorama_y,length_panorama_x,3));
x(1,1)=abs(x(1,1));
y(1,1)=abs(y(1,1));

if xini<=0
    if yini<=0
        img_panorama(1:yT,1:xT,:)=imTransf(:,:,:);
        img_panorama(y(1,1)+2:yB+y(1,1)+1,x(1,1)+2:xB+x(1,1)+1,:)=IBase(:,:,:);
    else
        img_panorama(y(1,1):yT+y(1,1)-1,1:xT,:)=imTransf(:,:,:);
        img_panorama(1:yB,x(1,1)+2:xB+x(1,1)+1,:)=IBase(:,:,:);
    end
else
    if yini<=0
        img_panorama(1:yT,x(1,1):xT+x(1,1)-1,:)=imTransf(:,:,:);
        img_panorama(y(1,1)+2:yB+y(1,1)+1,1:xB,:)=IBase(:,:,:);
    else
        img_panorama(y(1,1):yT+y(1,1)-1,x(1,1):xT+x(1,1)-1,:)=imTransf(:,:,:);
        img_panorama(1:yB,1:xB,:)=IBase(:,:,:);
    end
end

Panorama = img_panorama;
imBase = IBase;
imTransfor = imTransf;

figure;
subplot(2,2,1);subimage(ITransf);axis off
title('Imagem a Transformar')
subplot(2,2,3);subimage(imTransfor);axis off
title('Imagem Transformada')
subplot(2,2,4);subimage(imBase);axis off
title('Imagem Base')
subplot(2,2,2);subimage(Panorama);axis off
title('Panorama')

figure
imshow(img_panorama)

end
