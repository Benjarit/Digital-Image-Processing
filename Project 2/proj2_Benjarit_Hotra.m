clc; clear; close all; format compact;


input_image_file = 'Fig0333(a)(test_pattern_blurring_orig).tif'
image = imread(input_image_file);

figure(1);
subplot(2,2,1);
imshow(image);
caption = sprintf('Input Image');
title(caption, 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% section 1
avgFilterOutput = image;

% for 3x3 window
% N = 3;
% p = 1;

% for 5x5 window
N = 5;
p = 2;

const = 1/N^2;
avgFilter = const * ones(N,N);
inputImage1 = padarray(image,[p p],0, 'both');

[row,col,~] = size(inputImage1);

% Convolve Average filter with the original image
for i = 1+p:row-2*p
    for j = 1+p:col-2*p
        % element wise operation
        mul = avgFilter .* double(inputImage1(i-p:i+p,j-p:j+p));
        avgFilterOutput(i,j) = sum( sum( mul ));
    end
end

subplot(2,2,2);
imshow(avgFilterOutput);
caption = sprintf('Average Filter');
title(caption, 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussian fillter
GaussianFilteringOutput = image;

% for 3x3 window
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

subplot(2,2,3);
imshow(GaussianFilteringOutput);
caption = sprintf('Gussian Filter');
title(caption, 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Median filter
MedianFilteringOutput = image;

% for 3x3 window
N = 3;
p = 1;

% for 5x5 window
% N = 5;
% p = 2;

sigma = 1;

const = 1/(2*pi);
deno = -1/(2*(sigma^2));
middle = round(N/2);

sortedData = zeros(N,N);
medianFilter = zeros(N,N);
inputImage1 = padarray(image,[p p],0, 'both');

[row,col,~] = size(inputImage1);

% Convolve Median filter with the original image
for i = 1+p:row-2*p
    for j = 1+p:col-2*p
        % Rank its NxN neighbor
        medianInd = round(N*N/2);
        sortedData = sort(double(inputImage1(i-p:i+p,j-p:j+p)));
        
        medianNumber = sortedData(medianInd);
        MedianFilteringOutput(i,j) = medianNumber;
        
%         medianFilter(:,:) = medianNumber;
%         % element wise operation
%         mul = medianFilter .* double(inputImage1(i-p:i+p,j-p:j+p));
%         MedianFilteringOutput(i,j) = sum( sum( mul ));
    end
end

subplot(2,2,4);
imshow(MedianFilteringOutput);
caption = sprintf('Median Filter');
title(caption, 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sharpen 
input_image_file = 'Fig0338(a)(blurry_moon).tif';
image = imread(input_image_file);

figure(2);
subplot(1,2,1);
imshow(image);
caption = sprintf('Input Image');
title(caption, 'FontSize', 14);

SharpenImageOutput = image;

% for 3x3 window
% N = 3;
p = 1;

% var1 = [0 -1 0; -1 4 -1; 0 -1 0];
var2 = [-1 -1 -1; -1 8 -1; -1 -1 -1];
init = [0 0 0; 0 1 0; 0 0 0];
comb = var2 + init;


inputImage1 = padarray(image,[p p],0, 'both');

[row,col,~] = size(inputImage1);

% Convolve Median filter with the original image
for i = 1+p:row-2*p
    for j = 1+p:col-2*p        
        % element wise operation
        mul = comb .* double(inputImage1(i-p:i+p,j-p:j+p));
        SharpenImageOutput(i,j) = sum( sum( mul ));
    end
end
subplot(1,2,2);
imshow(SharpenImageOutput);
caption = sprintf('Sharpen Image');
title(caption, 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pyramid Gaussian

input_image_file = 'Lenna.png';
image = imread(input_image_file);

% convert to grayscale
image = rgb2gray(image);

sizeImage = [512 256 128 64 32 16 8];
currentImage = cell(1,size(sizeImage,2));

currentImage{1} = image;
downsampleImage = image;
L = size(sizeImage,2);

for k = 1:L
    
    % Gaussian filter
    % for 3x3 window
    N = 3;
    p = 1;
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

    inputImage1 = padarray(downsampleImage,[p p],0, 'both');

    [row,col,~] = size(inputImage1);

    % Convolve Gussian filter with the original image
    for i = 1+p:row-2*p
        for j = 1+p:col-2*p
            % element wise operation
            mul = gaussianFilter .* double(inputImage1(i-p:i+p,j-p:j+p));
            GaussianFilteringOutput(i,j) = sum( sum( mul ));
        end
    end % end of Gaussian Filter

    currentImage{k} = GaussianFilteringOutput; % store the downsample image
    downsampleImage = GaussianFilteringOutput(1:2:sizeImage(k), 1:2:sizeImage(k));
    GaussianFilteringOutput = downsampleImage;
end


figure(3);
subplot(1,7,1); 
imshow(currentImage{1});
caption = sprintf('512');
title(caption, 'FontSize', 14);

subplot(1,7,2); 
imshow(currentImage{2});
caption = sprintf('256');
title(caption, 'FontSize', 14);

subplot(1,7,3); 
imshow(currentImage{3});
caption = sprintf('128');
title(caption, 'FontSize', 14);

subplot(1,7,4); 
imshow(currentImage{4});
caption = sprintf('64');
title(caption, 'FontSize', 14);

subplot(1,7,5); 
imshow(currentImage{5});
caption = sprintf('32');
title(caption, 'FontSize', 14);

subplot(1,7,6); 
imshow(currentImage{6});
caption = sprintf('16');
title(caption, 'FontSize', 14);

subplot(1,7,7); 
imshow(currentImage{7});
caption = sprintf('8');
title(caption, 'FontSize', 14);