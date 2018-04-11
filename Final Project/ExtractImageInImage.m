clearvars; clc; format compact; close all;
im=imread('out.png');
figure();
imshow(im)

imsize  = size(im);
wtimage = zeros(imsize);

for i=1:imsize(1)
    for j=1:imsize(2)
        wtpixel = bitand(im(i,j), 1);
        if size(im,3) == 3 % color image
            if wtpixel == 0
                wtimage(i,j,1) = 0;
                wtimage(i,j,2) = 0;
                wtimage(i,j,3) = 0;
            else
                wtimage(i,j,1) = 255;
                wtimage(i,j,2) = 255;
                wtimage(i,j,3) = 255;
            end
        else % gray scale image
            if wtpixel == 0
                wtimage(i,j) = 0;
            else
                wtimage(i,j) = 255;
            end
        end
    end
end
figure();
imshow(wtimage)