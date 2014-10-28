[image1, descriptors1, locs1] = sift('img/rsz_imga.png');
[image2, descriptors2, locs2] = sift('img/rsz_imgb.png');

showkeys(image1, locs1);
showkeys(image2, locs2);

