clc; clear; close all; format compact;


input_image_file = 'Fig0309(a)(washed_out_aerial_image).tif'
image = imread(input_image_file);

figure(1);
subplot(1,2,1);
imshow(image);
caption = sprintf('Input Image');
title(caption, 'FontSize', 14);

% calculate the histogram of the input image
[row,col,num_channel] = size(image);

images_histogram = zeros(1,256);

%% calculate the histogram of the input image
for i = 1:256
    if (num_channel == 1)
        images_histogram(i) = sum(sum(image(:,:,1)==i-1));
    end

end 

% plot histrogram
figure(2);
subplot(2,1,1);
bar(images_histogram);
ylabel('frequency');
xlabel('pixel values');
title('Histogram');

%% perform cumulative histrogram
cumulative = zeros(1,256);

% initilaze the first position
cumulative(1) = images_histogram(1);

for i = 2:256
    cumulative(i) = cumulative(i-1) + images_histogram(i);
end

% plot cumulative histrogram
subplot(2,1,2);
bar(cumulative,'r');
ylabel('frequency');
xlabel('pixel values');
title('Cumulative Histogram');

%% perform histrogram equalization
I = image(:,:,1);
new_image = image(:,:,1);

if (num_channel == 1)

   % Numbers of pixels
   dimension = col * row;
  
   % Normalize cumulative histogram
   for i = 1:256
       
       % find() return matrix of pixel positions
       I = find(image(:,:,1)==i-1);
       
       % assign new value into new_image matriz
       new_image(I) = round(255*cumulative(i)/(dimension));
   end

   figure(1);
   subplot(1,2,2);
   imshow(new_image);
   caption = sprintf('Histogram Equalized Image');
   title(caption, 'FontSize', 14);
   
   % calculate the histogram of the new image
    for i = 1:256
       images_histogram(i) = sum(sum(new_image(:,:,1)==i-1));
    end
  
    % plot histrogram
    figure(3);
    subplot(2,1,1);
    bar(images_histogram);
    ylabel('frequency');
    xlabel('pixel values');
    title('Histogram');
    
    %% perform cumulative histrogram
    
    % initilaze the first position
    cumulative(1) = images_histogram(1);

    for i = 2:256
        cumulative(i) = cumulative(i-1) + images_histogram(i);
    end

    % plot cumulative histrogram
    subplot(2,1,2);
    bar(cumulative,'g');
    ylabel('frequency');
    xlabel('pixel values');
    title('Cumulative Histogram');

end