%% Video assignment, EE5120: Applied Linear Algebra I for EE

clear;
close all;
clc;

num_block = 6; % number of blocks to divide
num = 16; % number of top singular values to consider for each block
noise_level = 0.04; % set the noise mean for gaussian noise


Im = (imread('6.jpg'));
% imshow(Im);
% title('Original Rectangular Image');

Im = im2double(Im);
Im_gray = rgb2gray(Im); 
% imshow(Im_gray);
% title('Original Gray Scale Image');
%%

[rows, columns] = size(Im_gray); 
M = Im_gray' * Im_gray ;
[eig_vec,eigen] = eig(M); 
sigma = sqrt(eigen); 
C = Im_gray * eig_vec; 
U = zeros(rows,columns);
for i = 1:columns
    U(:,i) = C(:,i) / sigma(i,i);
end

Im_gray = imnoise(Im_gray,'gaussian',noise_level);
figure;
imshow(uint8(Im_gray*255))
title('Noised Image');

%% Break the image into blocks 

row_size = floor(rows/num_block);
col_size = floor(columns/num_block);

block = zeros(num_block^2,row_size,col_size);
a = 1;
for i = 0:num_block-1
    for l = 0:num_block-1
        for j = 1:row_size
            for k = 1:col_size
                block(a,j,k) = Im_gray(i *row_size + j , l*col_size + k);
            end
        end
        a = a + 1;
    end
end

% Display Sample block
% for j = 1:row_size
%     for k = 1:col_size
%         m(j,k) = block(36,j,k);
%     end
% end
% figure;
% imshow(uint8(m.*255))
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
for i = 0:num_block-1
    for l = 0:num_block-1
        for j = 1:row_size
            for k = 1:col_size
                Im_rec(i *row_size + j , l*col_size + k) = Im_hat(a,j,k);
            end
        end
        a = a + 1;
    end
end

figure; imshow(Im_rec,[])
title('Reconstructed Image')
