clc; clearvars; close all; format compact;

% This program is to implement the Canny edge detection algorithm
%% Filtered gradient
% 1. Load an image
input_image_file = 'Lenna.png';
image = imread(input_image_file);
image = rgb2gray(image);

figure(1);
subplot(1,2,1);
imshow(image);
caption = sprintf('Input Image');
title(caption, 'FontSize', 14);

%-------------------------------------------------------------------------%

% 2. Convolve the image with a Gaussian
GaussianFilteringOutput = image;

% % for 3x3 window
% N = 3;
% p = 1;

% for 5x5 window
N = 5;
p = 2;

sigma = 1;

const = 1/(2*pi);
deno = -1/(2*(sigma^2));
middle = round(N/2);

gaussianFilter = zeros(N,N);

for i = 1:N
    for j = 1:N
        gaussianFilter(i,j) = const * exp(deno * ((i-middle)^2 + (j-middle)^2));
    end
end

inputImage1 = padarray(image,[p p],0, 'both');

[row,col,~] = size(inputImage1);

% Convolve Gussian filter with the original image
for i = 1+p:row-2*p
    for j = 1+p:col-2*p
        % element wise operation
        mul = gaussianFilter .* double(inputImage1(i-p:i+p,j-p:j+p));
        GaussianFilteringOutput(i,j) = sum( sum( mul ));
    end
end

subplot(1,2,2);
imshow(GaussianFilteringOutput);
caption = sprintf('Gussian Filtering');
title(caption, 'FontSize', 14);

%-------------------------------------------------------------------------%
% 3. Find the x and y components of the gredient Fx and Fy at each point
Fx = image;
Fy = image;

% for 3x3 window
p = 1;
Gx = [-1 0 1; -2 0 2; -1 0 1];
Gy = Gx';


GaussianImage = padarray(GaussianFilteringOutput,[p p],0, 'both');
[row,col,~] = size(GaussianImage);

for i = 1+p:row-2*p
    for j = 1+p:col-2*p
        % element wise operation
        mul = Gx .* double(GaussianImage(i-p:i+p,j-p:j+p));
        Fx(i,j) = sum( sum( mul ));
        
        mu2 = Gy .* double(GaussianImage(i-p:i+p,j-p:j+p));
        Fy(i,j) = sum( sum( mu2 ));
    end
end

figure(2);
subplot(2,2,1);
imshow(Fx);
caption = sprintf('Fx');
title(caption, 'FontSize', 14);

subplot(2,2,2);
imshow(Fy);
caption = sprintf('Fy');
title(caption, 'FontSize', 14);

%-------------------------------------------------------------------------%
% 4. Compute the edge strength F and edge orientation D
Fx1 = double(Fx);
Fy2 = double(Fy);

F = hypot(double(Fx),double(Fy));
D = atan2(Fy2,Fx1);

subplot(2,2,3);
imshow(F);
caption = sprintf('F');
title(caption, 'FontSize', 14);

subplot(2,2,4);
imshow(D);
caption = sprintf('D');
title(caption, 'FontSize', 14);

%--------------------------------------------------
% Nonmaximum suppresion
% Create thinned edge image I(x,y)
D = (D * 180) / pi;
outOfMatrix = Inf;

oImg = F;
[r,c] = size(F);

for i=1:r %vertical
    for j=1:c %horizontal
        switch D(i,j)
            case 0
                %   0 - east
                east = outOfMatrix;
                west = outOfMatrix;
                if (j+1) <= c
                    east = F(i,j+1);
                end
                if (j-1) > 0
                    west = F(i,j-1);
                end 
                if F(i,j) < east || F(i,j) < west
                    oImg(i,j) = 0;
                end             
            case 1
                %   1 - north-east
                north_east = outOfMatrix;
                south_west = outOfMatrix;
                if (j+1) <= c && (i-1) > 0
                    north_east = F(i-1,j+1);
                end
                if (j-1) > 0 && (i+1) <= r
                    south_west = F(i+1,j-1);
                end 
                if F(i,j) < north_east || F(i,j) < south_west
                    oImg(i,j) = 0;
                end               
            case 2
                %   2 - north
                north = outOfMatrix;
                south = outOfMatrix;
                if (i-1) > 0
                    north = F(i-1,j);
                end
                if (i+1) <= r
                    south = F(i+1,j);
                end 
                if F(i,j) < north || F(i,j) < south
                    oImg(i,j) = 0;
                end                
            otherwise
                %   3 - north-west
                north_west = outOfMatrix;
                south_east = outOfMatrix;
                if (j-1) > 0 && (i-1) > 0
                    north_west = F(i-1,j-1);
                end
                if (j+1) <= c && (i+1) <= r
                    south_east = F(i+1,j+1);
                end 
                if F(i,j) < north_west || F(i,j) < south_east
                    oImg(i,j) = 0;
                end
        end
    end
end

oImg = oImg./max(max(oImg));

figure(3);
imshow(oImg);
caption = sprintf('Non-Maximum suppression');
title(caption, 'FontSize', 14);

% Hysteresis thresholding--------------------------------------------------
I_max = 30;
I_min = 20;
hys = hysthresh(oImg,I_max,I_min);
figure(4);
imshow(hys);
caption = sprintf('Hysteresis thresholding ');
title(caption, 'FontSize', 14);
% -----------------------------------------
function bw = hysthresh(im, I_max, I_min)
    if I_max < I_min
	tmp = I_max;
	I_max = I_min; 
	I_min = tmp;
    end
    aboveT2 = im > I_min;               
    [aboveT1r, aboveT1c] = find(im > I_max);  
    bw = bwselect(aboveT2, aboveT1c, aboveT1r, 8);
end