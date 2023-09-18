%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by ZhirongTang
% Date:     2022-05-06
% Default configurations
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
startup;

rf_freq = 470e6;
sf = 10;
bw = 125e3;
fs = 1e6;
cr = 4;
payload_len = 12;

% LoRaPHY
phy = LoRaPHY(rf_freq, sf, bw, fs, cr, payload_len);
phy.is_lowpass = true;
symbol_len = phy.symbol_len;
signal_len = phy.calc_samp_num(payload_len);    % 2*bw baseband sample points