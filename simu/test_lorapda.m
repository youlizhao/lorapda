%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by ZhirongTang
% Date:     2022-05-08
% Revised:  2023-08-30
% Func: Comparing the SNR-BER of standard LoRa decoder(hard-decision decoding) 
%       and LoRaPPA soft decoder(soft-decision decoding)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seed = 10000;
rng(seed);

% phy configs
real_config;

% evaluation params
num_user = 2;
max_cfo = 5e3;
max_to = round(0.1*2^sf);
signal_len = signal_len + max_to*2;
snr = 10;
iter_len = 10;

% result
num_user_range = num_user;
snr_range = snr;
ser_ppa = -ones(iter_len, length(snr_range), length(num_user_range));
ber_hard = -ones(size(ser_ppa));
ber_soft = -ones(size(ser_ppa));
sum_ber_hard = -ones(size(ser_ppa));
sum_ber_soft = -ones(size(ser_ppa));

% topk probability for LoRaPPA-soft
topk = 2;

tic;
% LoRaPPA
cut_to = max_to + 2;
threshold = 5;
fast_mode = false;
ppa = LoRaPPA(rf_freq, sf, bw, fs, cr, payload_len, num_user, cut_to, threshold, fast_mode);
ppa.topk = topk;
ppa.is_real = false;
ppa.is_resample = false;

% phy params
amp_range = [0.6; 0.9; 1.3; 1.7; 2; 2.3];
amp = amp_range(1: num_user);

for iter = 1: iter_len   % may use PARFOR here
    %%%%%%%% init: gen cfo/to
    to = unifrnd(0, max_to, [num_user, 1]);    % to with frac
    phase = exp(1j*2*pi*rand(num_user, 1));
    
    % groud truth (symbols&payload, sfo_drift, amplitude)
    payload_tx = randi([0, 255], payload_len, num_user);
    symbols_tx = zeros(symbol_len, num_user);
    cfo = unifrnd(-max_cfo, max_cfo, [num_user, 1]);
    cfo_shift = cfo*2^sf/bw;
    for ui = 1: num_user
        symbols_tx(:, ui) = phy.encode(payload_tx(:,ui));

    end
    %%%%%%%%
    
    % add to & noise
    signal_tx = zeros(signal_len, 1);
    for ui = 1: num_user
        sig = phy.real_modulate(symbols_tx(:, ui), amp(ui), phase(ui), cfo(ui), to(ui));
        sig = resample(sig, 2*bw, fs);
        signal_tx(1: length(sig)) = signal_tx(1: length(sig)) + sig;
    end
    % add noise
    signal_rx = utils.add_noise(signal_tx, snr);
    
    % demodulation
    ppa.preamble_est(signal_rx);
    [symbols_ppa, symbols_ppa_p0] = ppa.demodulate();
    
    % decode
    payload_hard = zeros(payload_len, num_user);
    payload_soft = zeros(size(payload_hard));
    for ui = 1: num_user
        payload_hard(:, ui) = phy.decode(symbols_ppa(:, ui));
        payload_soft(:, ui) = ppa.soft_decode(symbols_ppa_p0(:, ui));
    end
    
    % ser, ber calculation
    snr_i = 1;
    num_user_i = 1;
    ser_ppa(iter, snr_i, num_user_i) = utils.calc_ser(symbols_ppa, symbols_tx);
    ber_hard(iter, snr_i, num_user_i) = utils.calc_ber(payload_hard, payload_tx);
    ber_soft(iter, snr_i, num_user_i) = utils.calc_ber(payload_soft, payload_tx);
    sum_ber_hard(iter, snr_i, num_user_i) = utils.calc_sum_ber(payload_hard, payload_tx);
    sum_ber_soft(iter, snr_i, num_user_i) = utils.calc_sum_ber(payload_soft, payload_tx);
    
    fprintf("\niter:%d, snr:%d, num_user:%d\n", iter, snr, num_user);
    fprintf("ser_ppa:%g\n", ser_ppa(iter, snr_i, num_user_i));
    fprintf("ber_hard:%g, ber_soft:%g\n", ber_hard(iter, snr_i, num_user_i), ber_soft(iter, snr_i, num_user_i));
    fprintf("sum_ber_hard:%g, sum_ber_soft:%g\n", sum_ber_hard(iter, snr_i, num_user_i), sum_ber_soft(iter, snr_i, num_user_i));
end     % end of iter
toc;