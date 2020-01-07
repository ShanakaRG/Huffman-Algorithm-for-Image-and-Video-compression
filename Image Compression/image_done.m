close all;
clear all;
clc;

OriginalImage='sample_1.jpg';
img_rgb=imread(OriginalImage); % read the original image
img = rgb2gray(img_rgb); %convert original image to Gray scale
R_plane=img_rgb(:,:,1);
G_plane=img_rgb(:,:,2);
B_plane=img_rgb(:,:,3);
back_to_original = cat(3,R_plane,G_plane,B_plane);
figure;imshow(img_rgb);title('orininal Image');
figure;imshow(R_plane);title('Red plane Image');
figure;imshow(G_plane);title('Green plane Image');
figure;imshow(B_plane);title('Blue plane Image');
figure;imshow(back_to_original);title('Back to original Image');

figure;imshow(img);title('GRAY Image');
[r c b]=size(img);



 
fun = @dct2;
dct_gray = blkproc(img,[8 8],fun);
dct_red = blkproc(R_plane,[8 8],fun);
dct_green = blkproc(G_plane,[8 8],fun);
dct_blue= blkproc(B_plane,[8 8],fun);
figure;imshow(dct_gray);title('DCT Image GRAY');

fun1 = @fact;
quantized_gray = blkproc(dct_gray,[8 8],fun1);
quantized_red = blkproc(dct_red,[8 8],fun1);
quantized_green = blkproc(dct_green,[8 8],fun1);
quantized_blue = blkproc(dct_blue,[8 8],fun1);
figure;imshow(quantized_gray);title('quantized Image GRAY')


round_mat_gray = round(quantized_gray);
round_mat_red = round(quantized_red);
round_mat_green = round(quantized_green);
round_mat_blue = round(quantized_blue);

transpose_gray=round_mat_gray';
transpose_red=round_mat_red';
transpose_green=round_mat_green';
transpose_blue=round_mat_blue';

column_gray=transpose_gray(:);
column_red=transpose_red(:);
column_green=transpose_green(:);
column_blue=transpose_blue(:);

[g_gray,~,Intensity_val_gray] = grp2idx(column_gray);
Frequency_gray  = accumarray(g_gray,1);

[g_red,~,Intensity_val_red] = grp2idx(column_red);
Frequency_red  = accumarray(g_red,1);

[g_green,~,Intensity_val_green] = grp2idx(column_green);
Frequency_green  = accumarray(g_green,1);

[g_blue,~,Intensity_val_blue] = grp2idx(column_blue);
Frequency_blue  = accumarray(g_blue,1);



[Intensity_val_gray Frequency_gray];
probability_gray=Frequency_gray./sum(Frequency_gray); % probabilities

[Intensity_val_red Frequency_red];
probability_red=Frequency_red./sum(Frequency_red); % probabilities

[Intensity_val_green Frequency_green];
probability_green=Frequency_green./sum(Frequency_green); % probabilities

[Intensity_val_blue Frequency_blue];
probability_blue=Frequency_blue./sum(Frequency_blue); % probabilities

% T = table(Intensity_val_gray,Frequency_gray,probability_gray); %table (element | count | probability )
% T(1:length(Frequency_gray),:);

DICT_gray = huffmandict(Intensity_val_gray, probability_gray) ;
DICT_red = huffmandict(Intensity_val_red, probability_red) ;
DICT_green = huffmandict(Intensity_val_green, probability_green) ;
DICT_blue = huffmandict(Intensity_val_blue, probability_blue) ;

huff_encoded_gray = huffmanenco(column_gray,DICT_gray);
huff_encoded_red = huffmanenco(column_red,DICT_red);
huff_encoded_green = huffmanenco(column_green,DICT_green);
huff_encoded_blue = huffmanenco(column_blue,DICT_blue);

total_size = size(huff_encoded_red)+size(huff_encoded_green)+size(huff_encoded_blue)
original_image = 640*400*3*8
ratio = original_image/ total_size(1)

huff_deco_gray = huffmandeco(huff_encoded_gray,DICT_gray);
huff_deco_red = huffmandeco(huff_encoded_red,DICT_red);
huff_deco_green = huffmandeco(huff_encoded_green,DICT_green);
huff_deco_blue = huffmandeco(huff_encoded_blue,DICT_blue);

% equal_1 = isequal(quant8_I(:),deco) % Check whether the decoding is correct.
deco_img_gray = vec2mat(huff_deco_gray,c);
deco_img_red = vec2mat(huff_deco_red,c);
deco_img_green = vec2mat(huff_deco_green,c);
deco_img_blue = vec2mat(huff_deco_blue,c);

figure;imshow (deco_img_gray,[]);title('Decoded Image GRAY')
figure;imshow (deco_img_red,[]);title('Decoded Image red')
figure;imshow (deco_img_green,[]);title('Decoded Image green')
figure;imshow (deco_img_blue,[]);title('Decoded Image blue')
equal_2 = isequal(round_mat_gray,deco_img_gray)


fun2 = @fact2;
dequantized_gray = blkproc(deco_img_gray,[8 8],fun2);
dequantized_red = blkproc(deco_img_red,[8 8],fun2);
dequantized_green = blkproc(deco_img_green,[8 8],fun2);
dequantized_blue = blkproc(deco_img_blue,[8 8],fun2);

figure;imshow(dequantized_gray);title('dequantized Image GRAY')
figure;imshow(dequantized_red);title('dequantized Image red')
figure;imshow(dequantized_green);title('dequantized Image green')
figure;imshow(dequantized_blue);title('dequantized Image blue')

fun22=@idct2;
inv_DCT_gray = blkproc(dequantized_gray , [8 8],fun22);
inv_DCT_red= blkproc(dequantized_red , [8 8],fun22);
inv_DCT_green = blkproc(dequantized_green , [8 8],fun22);
inv_DCT_blue = blkproc(dequantized_blue , [8 8],fun22);

reconstruct_original_img = cat(3, inv_DCT_red, inv_DCT_green, inv_DCT_blue);
img_recieved = uint8(reconstruct_original_img);
figure;imshow(inv_DCT_gray,[]); title('inverse DCT Image GRAY')
figure;imshow(img_recieved); title('inverse DCT Image RGB / Recieved image')




function f = fact(n)
    f = (1/2^1).*n;
    
end

function f2 = fact2(m)
    f2 = 2^1*m;
end








