clearvars; clc; format compact; close all;
% wtmark function performs watermarking
% im     = Input Image
% wt     = Watermark
% embimg = Output Embedded image

% Read in image
im = imread('b.jpg');

% Read in watermark
wt = imread('w.jpg');


figure();
imshowpair(im, wt, 'montage');

%% Resize watermark (pad if smaller and resize if bigger)
imsize = size(im); % get size of image
wtsize = size(wt); % get size of watermark

% Resizing row 
scaleX = imsize(1)/wtsize(1);
if imsize(1) < wtsize(1)
    wt = imresize(wt, [scaleX NaN]);
end

% Resizing column
scaleY = imsize(2)/wtsize(2);
if imsize(2) < wtsize(2)
    wt = imresize(wt, [scaleY NaN]);
end

wt = imresize(wt, [size(im,1) NaN]);
 
figure();
imshowpair(im, wt, 'montage')

%% RGB 
for i=1:imsize(1) 
    for j=1:imsize(2)
        if size(im,3) == 3
            wtpixel = zeros(3,1);
            wtpixel(1) = wt(i,j,1);
            wtpixel(2) = wt(i,j,2);
            wtpixel(3) = wt(i,j,3);
            if wtpixel >= 128
                im(i,j,1) = bitor(im(i,j,1), 1);
                im(i,j,2) = bitor(im(i,j,2), 1);
                im(i,j,3) = bitor(im(i,j,3), 1);
            else
                im(i,j,1) = bitand(im(i,j,1), 254);
                im(i,j,2) = bitand(im(i,j,2), 254);
                im(i,j,3) = bitand(im(i,j,3), 254);
            end
        else
            wtpixel = wt(i,j);
            if wtpixel >= 128
                im(i,j) = bitor(im(i,j), 1);
            else
                im(i,j) = bitand(im(i,j), 254);
            end
        end
    end
end
figure();
imshow(im)
imwrite(im, 'out.png');