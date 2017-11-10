%% Video assignment, EE5120: Applied Linear Algebra I for EE
% Author: Honey Gupta, MS Scholar
% EE Department, IIT Madras


clear;
close all;
clc;

num_block = 6; % number of blocks to divide
num = 16; % number of top singular values to consider for each block
noise_level = 0.04; % set the noise mean for gaussian noise


Im = (imread('rect_6.jpg'));
imshow(Im);
title('Original Rectangular Image');

Im = im2double(Im);
Im_gray = rgb2gray(Im); 
figure;
imshow(Im_gray);
title('Original Gray Scale Image');
%% Calculating SVD and adding noise to 

Im_gray = imnoise(Im_gray,'gaussian',noise_level);
figure;
imshow(uint8(Im_gray*255))
title('Noised Image (Added gaussian noise with mean = 0.04, variance =  0.01)');

[rows, columns] = size(Im_gray); 
M = Im_gray' * Im_gray ;
[eig_vec,eigen] = eig(M); 
sigma = sqrt(eigen); 
C = Im_gray * eig_vec; 
U = zeros(rows,columns);
for i = 1:columns
    U(:,i) = C(:,i) / sigma(i,i);
end


%% Break the image into blocks 

row_size = floor(rows/num_block);
col_size = floor(columns/num_block);

block = zeros(num_block^2,row_size,col_size);
a = 1;
figure;
title('Segments of Noisy image');
for i = 0:num_block-1
    for l = 0:num_block-1
        for j = 1:row_size
            for k = 1:col_size
                block(a,j,k) = Im_gray(i *row_size + j , l*col_size + k);
                m(j,k) = Im_gray(i *row_size + j , l*col_size + k);

            end
        end
        subplot(num_block,num_block,a)
        imshow(uint8(m.*255))
        a = a + 1;
    end
end

%% Perform SVD on the blocks

Im_hat = zeros(size(block));
for i = 1:num_block*num_block
    for j = 1:row_size
        for k = 1:col_size
            m(j,k) = block(i,j,k);
        end
    end
    M = m' * m ;
    [eig_vec,eigen] = eig(M); 
    sigma = sqrt(eigen); 
    C = m * eig_vec; 
    U = zeros(row_size,col_size);
    for j = 1:col_size
        U(:,j) = C(:,j) / sigma(j,j);
    end
    
    
    sigma_re = zeros(col_size,col_size);
    for j=size(sigma,2) - num : size(sigma,2)
        sigma_re(j,j) = sigma(j,j);
    end

    Im_hat(i,:,:) =( U * sigma_re * eig_vec') .* 255; %Reconstructed image
    
end


a = 1;
figure;
title('Segments of Reconstructed Image');
clear m;
for i = 0:num_block-1
    for l = 0:num_block-1
        for j = 1:row_size
            for k = 1:col_size
                Im_rec(i *row_size + j , l*col_size + k) = Im_hat(a,j,k);
                m(j,k) = Im_hat(a,j,k);
            end
        end
        subplot(num_block,num_block,a)
        imshow(uint8(m))
        a = a + 1;
    end
end


figure; 
subplot(2,1,1)
imshow(Im_gray,[])
title('Original Image');
subplot(2,1,2)
imshow(Im_rec + 20,[]) % added 20 to compensate for the intensity loss
title('Reconstructed Image')
