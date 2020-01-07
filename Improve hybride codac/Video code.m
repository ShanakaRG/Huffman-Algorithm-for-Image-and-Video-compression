clear all;
close all;
clc;
no_of_frames=50;
Quality = 0.1
p=7
mbSize=8
t=5
total_size = [0 0]
original_size = 0
[Y U V]= yuv_import('soccer_cif.yuv',[352 288],no_of_frames);
%initialize the motion vector
k=1;

C=xcorr2(Y{1,1},Y{1,2});
figure ; mesh (C);title ('Correlation')
figure; plot(C(:));

for j =1:t:no_of_frames-1
motion_prediction_Y = zeros(288,352);
motion_prediction_u = zeros(144,176);
motion_prediction_v = zeros(144,176);

predict_out_Y = zeros(288,352);
predict_out_U = zeros(144,176);
predict_out_V = zeros(144,176);

for i=1:1:t

frm_y = Y{1,k};
frm_u = U{1,k};
frm_v = V{1,k};



error_y = frm_y-motion_prediction_Y; 
error_u = frm_u-motion_prediction_u;
error_v = frm_v-motion_prediction_v;

[r, c, b]=size(frm_y);
[r_u, c_u, b_u]=size(frm_u);

%apply DCT for 8*8 blocks inthe image
fun = @dct2;
dct_y = blkproc(error_y,[8 8],fun);
dct_u = blkproc(error_u,[8 8],fun);
dct_v = blkproc(error_v,[8 8],fun);

% figure;imshow(dct_y);title('DCT Image_GRAY');
q_mtx = [16 11 10 16 24 40 51 61; 
         12 12 14 19 26 58 60 55;
         14 13 16 24 40 57 69 56; 
         14 17 22 29 51 87 80 62;
         18 22 37 56 68 109 103 77;
         24 35 55 64 81 104 113 92;
         49 64 78 87 103 121 120 101;
         72 92 95 98 112 100 103 99]*Quality;
%quantize the image

quantized_y = blkproc(dct_y,[8 8],'round(x./P1)',q_mtx);
quantized_u = blkproc(dct_u,[8 8],'round(x./P1)',q_mtx);
quantized_v = blkproc(dct_v,[8 8],'round(x./P1)',q_mtx);

% figure;imshow(quantized_y);title('Quantized Image GRAY')
% round_mat_y{1,i} = round(quantized_y); % round the quantized matrix
% round_mat_u{1,i} = round(quantized_u);
% round_mat_v{1,i} = round(quantized_v);

%dequatized and invDCT to get decode image

dequntize_y = blkproc(quantized_y,[8 8],'x.*P1',q_mtx);
dequntize_u = blkproc(quantized_u,[8 8],'x.*P1',q_mtx);
dequntize_v = blkproc(quantized_v,[8 8],'x.*P1',q_mtx);

fun = @idct2;
INVdct_y = blkproc(dequntize_y,[8 8],fun);
INVdct_u = blkproc(dequntize_u,[8 8],fun);
INVdct_v = blkproc(dequntize_v,[8 8],fun);

store_y{1,i} = INVdct_y + motion_prediction_Y;
store_u{1,i} = INVdct_u + motion_prediction_u;
store_v{1,i} = INVdct_v + motion_prediction_v;

[MV_Y,  DScomputations_y] = motionEstDS(Y{1,i+1},store_y{1,i}, mbSize, p);
[MV_u,  DScomputations_u] = motionEstDS(U{1,i+1},store_u{1,i}, mbSize, p);
[MV_v,  DScomputations_v] = motionEstDS(V{1,i+1},store_v{1,i}, mbSize, p);

motion_prediction_Y = motionComp(store_y{1,i}, MV_Y, mbSize);
motion_prediction_u = motionComp(store_u{1,i}, MV_u, mbSize);
motion_prediction_v = motionComp(store_v{1,i}, MV_v, mbSize);

[DICT_y,huff_quantized_y]= Henco(quantized_y);
[DICT_u,huff_quantized_u]= Henco(quantized_u);
[DICT_v,huff_quantized_v]= Henco(quantized_v);

[DICT_MV_y,huff_mv_y]= Henco(MV_Y);
[DICT_MV_u,huff_mv_u]= Henco(MV_u);
[DICT_MV_v,huff_mv_v]= Henco(MV_v);
%++++++++++++++++++++++++++++++++++++++++ DECODER  +++++++++++++++++++++++++++++
% %dequantization of the decoded matrix

huff_deco_quantized_y = huffmandeco(huff_quantized_y,DICT_y); %haffman decoded
huff_deco_quantized_u = huffmandeco(huff_quantized_u,DICT_u);
huff_deco_quantized_v = huffmandeco(huff_quantized_v,DICT_v);

huff_deco_MV_y = huffmandeco(huff_mv_y,DICT_MV_y);
huff_deco_MV_u = huffmandeco(huff_mv_u,DICT_MV_u);
huff_deco_MV_v = huffmandeco(huff_mv_v,DICT_MV_v);

total_size = total_size + size(huff_quantized_y)+size(huff_quantized_u)+size(huff_quantized_v)+size(huff_mv_y)+size(huff_mv_u)+size(huff_mv_v)
original_size = original_size + ((288*352)+2*(144*176))*8

deco_qunt_y{1,k} = vec2mat(huff_deco_quantized_y,c); %convert decode vector to matrix
deco_qunt_u{1,k} = vec2mat(huff_deco_quantized_u,c/2);
deco_qunt_v{1,k} = vec2mat(huff_deco_quantized_v,c/2);

deco_mv_y{1,k} = vec2mat(huff_deco_MV_y,1584); %convert decode vector to matrix
deco_mv_u{1,k} = vec2mat(huff_deco_MV_u,396);
deco_mv_v{1,k} = vec2mat(huff_deco_MV_v,396);

    
Dequntize_y=blkproc(deco_qunt_y{1,k},[8 8],'x.*P1',q_mtx);  
Dequntize_u=blkproc(deco_qunt_u{1,k},[8 8],'x.*P1',q_mtx);  
Dequntize_v=blkproc(deco_qunt_v{1,k},[8 8],'x.*P1',q_mtx); 
    
fun1=@idct2;
Y_out{1,k}=blkproc(Dequntize_y,[8 8],fun1);
U_out{1,k}=blkproc(Dequntize_u,[8 8],fun1);
V_out{1,k}=blkproc(Dequntize_v,[8 8],fun1);

Final_Y{1,k}=Y_out{1,k}+predict_out_Y;
Final_U{1,k}=U_out{1,k}+predict_out_U;
Final_V{1,k}=V_out{1,k}+predict_out_V;
    
predict_out_Y=motionComp(Final_Y{1,k},deco_mv_y{1,k}, mbSize);
predict_out_U=motionComp(Final_U{1,k},deco_mv_u{1,k}, mbSize);
predict_out_V=motionComp(Final_V{1,k},deco_mv_v{1,k}, mbSize);

k=k+1
end
end
ratio = original_size/ total_size (1)
% yuv_export(rx_y,rx_u,rx_v,'akiyo_cif_output.yuv',no_of_frames)
yuv_export(Final_Y,Final_U,Final_V,'DECODE_OUTPUT_QP3.yuv',no_of_frames-1);
%function for quatization
function f = fact(n)
    f = (1/2^3).*n;
end
%function for dequatization
function f2 = fact2(m)
    f2 = (2^3).*m;
end

function [DICT, huff_encoded] = Henco(m)
transpose=m'; % get transpose
column=transpose(:);% convert matrix to a column vector
%probability calculation
[g,~,Intensity_val] = grp2idx(column);
Frequency  = accumarray(g,1);
[Intensity_val Frequency];
probability=Frequency./sum(Frequency); % probabilities
% buld  the huffman dictionary 
DICT = huffmandict(Intensity_val, probability) ; 
huff_encoded = huffmanenco(column,DICT); % huffman encoded
% huff_deco = huffmandeco(huff_encoded,DICT); %haffman decoded
end

function inv_DCT = iDCT(deco_img)
%dequantization of the decoded matrix
fun2 = @fact2;
inv_Quantize = blkproc(deco_img,[8 8],fun2);
Idct=@idct2;
inv_DCT = blkproc(inv_Quantize , [8 8],Idct);
end
