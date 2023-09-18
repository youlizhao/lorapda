%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by ZhirongTang
% Date:     2022-03-14
% Update:   2022-05-06
% LoRa Physical Layer Packet Aggregation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef LoRaPPA < handle & matlab.mixin.Copyable & LoRaPHY
    properties
        num_user                % number of users
        cut_to                  % # of baseband samples to be cut
        cfo_u                   % carrier frequency BIN offset of all users
        to_u                    % time BIN offset of all users
        chan_u                  % signal channel of all users
        sig_start               % preamble start indexes on baseband signal
        pream_est               % `coarse`: preamble coarse estimation
                                % `fine`: preamble fine estimation
        
        algo                    % demodulate algorithm
                                % `mlp`: maximum likelihood permutation
                                % `unc`: unconstrained optimization
                                % `bound`: full-search, as upper bound
        enum_algo               % enumeration algorithm
                                % `Match-M!`: peak matching first, then full peak enumeration for those who failed to match
                                % `M!`: 1~M full peak enumeration
                                % `M^M`: M peak enumeration
                                % `V^M`: V peak enumeration (N > V >= M)
        topk                    % keep topk possible permutation for sdd
        win                     % window type
        win_beta                % window with shape factor beta
        threshold               % peak higher than threshold is considered
        
        upchirp_corr            % upchirp for correlation
        upchirp_c               % upchirp after cut TO
        downchirp_c             % downchirp after cut TO
        upchirp_cw              % upchirp after cut TO and add window
        downchirp_cw            % downchirp after cut TO and add window
        fs_constr               % sampling frequency of reconstruction signal
                                % large fs leads to more precise and slower construction
        num_constr_samples      % number of sample points per reconstructed symbol
        fast_mode               % set `true` to use fast reconstruction mode
        fast_chirp              % pre-reconstructed base upchirps
        
        is_resample             % set `true` to resample signals to 2*bw
        is_real                 % set `true` for real signal process, which doesn't issue errors
        sfo_drift               % sampling frequency offset drift of real signals

        debug_ppa               % set `true` for debug information
        cfo_coarse_u            % coarse-grained carrier frequency BIN offset of all users
        to_coarse_u             % coarse-grained time BIN offset of all users
        chan_coarse_u           % coarse-grained signal channel of all users
    end
    
    methods
        function self = LoRaPPA(rx_freq, sf, bw, fs, cr, payload_len, num_user, cut_to, threshold, fast_mode)
            % LoRaPPA  Constructor
            self = self@LoRaPHY(rx_freq, sf, bw, fs, cr, payload_len);
            
            if nargin == 8
                self.threshold = 5;
                self.fast_mode = false;
            elseif nargin == 9
                self.threshold = threshold;
                self.fast_mode = false;
            elseif nargin == 10
                self.threshold = threshold;
                self.fast_mode = fast_mode;
            end
            
            self.num_user = num_user;
            self.cut_to = cut_to;
            self.enum_algo = "M^M";
            self.topk = 2;
            
            self.cfo_u = zeros(num_user, 1);
            self.to_u = zeros(num_user, 1);
            self.chan_u = zeros(num_user, 1);
            self.pream_est = "coarse";
            
            self.algo = "mlp";
            self.win_beta = 5;
            self.win = kaiser(self.num_samples, self.win_beta);
            
            self.debug_ppa = false;
            self.is_resample = true;
            self.is_real = false;
            
            self.sfo_drift = zeros(self.num_user, self.symbol_len);
            
            self.init_ppa();
        end
        
        function init_ppa(self)
            % init_ppa  Initialize some parameters
            self.fs_constr = 10*self.bw;
            self.num_constr_samples = 10*self.num_bb_samples;
            
            self.upchirp_corr = LoRaPPA.chirp(self.sf, self.bw, self.fs, 0, true);
            
            self.upchirp_c = self.upchirp;
            self.upchirp_c(1: 2*self.cut_to) = 0;
            self.upchirp_cw = self.upchirp_c .* self.win;
            self.downchirp_c = self.downchirp;
            self.downchirp_c(1: 2*self.cut_to) = 0;        
            self.downchirp_cw = self.downchirp_c .* self.win;
            
            if self.fast_mode
                self.fast_chirp = zeros(self.num_bb_samples*(self.fs_constr/self.bw), self.num_bb_samples);
                for si = 0: self.num_bb_samples - 1
                   self.fast_chirp(:, si+1) = LoRaPPA.chirp(self.sf, self.bw, self.fs_constr, si, true); 
                end
            end
        end
        
        function [symbols_m, payload_hard, payload_soft] = real_demodulate(self, raw_sig, algo)
            % real_demodulate  LoRa aggregated packet demodulation for real
            % signal
            %
            % input:
            %     raw_sig: raw signals
            %     algo: `mlp`, `unc`, `bound`
            % output:
            %     symbols_m: A matrix containing the demodulated results.
            %     payload_hard: A matrix containing decoded data with
            %                   standard lora hard decoder
            %     payload_soft: A matrix containing decoded data with
            %                   lorappa soft decoder
            
            if nargin == 3
                self.algo = algo;
            end
            self.is_real = true;
            
            % preamble detection
            self.detect(raw_sig);

            % channel & offset estimation
            self.preamble_est();

            % demodulation
            [symbols_m, symbols_p0] = self.demodulate();

            % soft decoding
            payload_hard = zeros(self.payload_len, self.num_user);
            payload_soft = zeros(self.payload_len, self.num_user);
            for ui = 1: self.num_user
                payload_hard(:, ui) = self.decode(symbols_m(:, ui));
                payload_soft(:, ui) = self.soft_decode(symbols_p0(:, ui));
            end
        end
        
        function [symbols_m, symbols_p0] = demodulate(self, algo, enum_algo, sig, varargin)
            % demodulate  LoRa aggregated packet demodulation
            %
            % input:
            %     algo: `mlp`, `unc`, `bound`
            % output:
            %     symbols_m: A matrix containing the demodulated results.
            %                Each column vector represents the symbols of
            %                a successfully demodulated packet.
            %     symbols_p0: A matrix containing probabilities of
            %                demodulated bits that equal to zero.
            %%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % global amp;
            % global phase;
            % global cfo;
            % global to; 
            % global symbols_tx;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % NOTE: perform preamble_est() first!       
            if nargin >= 2
                self.algo = algo;
            end
            if nargin >= 3
               self.enum_algo = enum_algo; 
            end
            if nargin >= 4
                self.sig = sig;
            end
            if nargin >= 5
                amp = varargin{1};
                phase = varargin{2};
                cfo = varargin{3};
                to = varargin{4};
                symbols_tx = varargin{5};
                self.cfo_u = cfo * self.num_bb_samples / self.bw;
                self.to_u = -abs(to);
                self.chan_u = amp.*phase;
            end
            
            % symbol demodulation
            sym_shift = self.cfo_u + self.to_u;
            
            symbols_m = zeros(self.symbol_len, self.num_user);
            symbols_p0 = zeros(self.symbol_len*self.sf, self.num_user);
            pl_start = round((self.preamble_len+4.25)*self.num_samples)+1;
            % rectangular window fft
            [pk, ft] = self.dechirp(pl_start, self.symbol_len, self.downchirp_c);
            pk(:, 2, :) = (pk(:, 2, :)) / self.zero_padding_ratio;          % FIXME: donnot need -1?
            % kaiser window fft
            [pk_w, ft_w] = self.dechirp(pl_start, self.symbol_len, self.downchirp_cw);
            pk_w(:, 2, :) = (pk_w(:, 2, :)) / self.zero_padding_ratio;      % FIXME: donnot need -1?
            
            % cfo start phase compenstation
            % NOTE: cfo start should be calculated on baseband samples
            % NOTE: cfo starts with 0
            cfo_start = (self.preamble_len+4.25)*self.num_bb_samples;
            
            for si = 1: self.symbol_len
                % remove -1 padding
                pk_ = pk_w(:,2,si);
                pk_(pk_<0) = [];
                
                % permutations of valid peaks
                sym_perm = self.peak_enum(pk_, self.num_user, self.enum_algo, pk_w(:, 1, si)).' - sym_shift;
                % use round() for fast_mode
                % FIXME: seems that fast_mode degrades the performance
                sym_perm = mod(sym_perm, self.num_bb_samples);
                
                if self.algo == "mlp" || self.algo == "bound"
                    % maximum likelihood permutation algo
                    perm_dist = zeros(size(sym_perm, 2), 1);
                    % may use PARFOR here
                    for pi = 1: size(sym_perm, 2)
                        % compute residual value for each permutation
                        sym = sym_perm(:, pi);
                        % perm_dist(pi) = LoRaPPA.residual_func(ft(:,si), self.sig_constr("freq", sym, self.cfo_u, self.to_u, true, cfo_start+(si-1)*self.num_bb_samples, self.chan_u), "freq");
                        % calculate residual function on top bins, which
                        % can improve soft-decoding performance
                        perm_dist(pi) = LoRaPPA.residual_func_bin(ft(:,si), self.sig_constr("freq", sym, self.cfo_u, self.to_u, true, cfo_start+(si-1)*self.num_bb_samples, self.chan_u), "freq", pk_*self.zero_padding_ratio);
                    end
                    
                    % hard-decision: keep minimum distance
                    [~, min_i] = min(perm_dist);
                    symbol = mod(round(sym_perm(:, min_i) - self.sfo_drift(:,si)), self.num_bb_samples);
                    
                    % soft-decision: keep topk minimum distance
                    [~, min_k] = mink(perm_dist, self.topk);
                    topk_perm = mod(round(sym_perm(:, min_k).' - self.sfo_drift(:,si).'), self.num_bb_samples);       % row: topk, col: user
                    topk_dist = perm_dist(min_k);
                    % convert per-symbol probability to per-bit probability
                    for ui = 1: self.num_user
                        % binary perm
                        perm_bin = de2bi(topk_perm(:, ui), self.sf, 2, 'left-msb');
                        per_bit_dist = zeros(2, self.sf);
                        % copy dist to per-bit
                        for bi = 1: self.sf
                           % bin == 0
                           per_bit_dist(1, bi) = sum(topk_dist(perm_bin(:, bi) == 0));
                           % bin == 1
                           per_bit_dist(2, bi) = sum(topk_dist(perm_bin(:, bi) == 1));
                        end
                        per_bit_dist(per_bit_dist == 0) = realmax();
                        % normalization
                        per_bit_dist = per_bit_dist ./ sum(per_bit_dist, 1);
                        % dist->pr
                        per_bit_dist = flipud(per_bit_dist);  
                        % keep the probability of 0
                        symbols_p0((si-1)*self.sf+1: si*self.sf, ui) = per_bit_dist(1, :);
                    end
                else
                    error('Error. Invalid algo.');
                end
                
                symbol = mod(round(symbol), self.num_bb_samples);
                symbols_m(si, :) = symbol;
                
            end
        end
        
        function preamble_est(self, sig)
            % preamble_est  Joint estimate cfo, to, channel of all users
            %
            % input:
            %     sig: input raw signal
            %

            if nargin == 2
                self.sig = sig;
            end
            
            % resample signal with 2*bandwidth
            if self.is_resample
                if self.is_lowpass
                    sig = lowpass(sig, self.bw/2, self.fs);
                end
                % self.sig = resample(sig, 2*self.bw, self.fs);
                self.sig = sig(1: self.fs/(2*self.bw): end);
            end
            
            % fft on both rectangular and kaiser window to find valid peaks
            % preamble upchirp shift
            [pk_u, ft_u] = self.dechirp(1, self.preamble_len, self.downchirp_cw);
            pk_u(:, 2, :) = (pk_u(:, 2, :) - 1) / self.zero_padding_ratio;
            if self.debug_ppa
                pk_w = self.dechirp(1, self.preamble_len, self.downchirp_cw);
                 pk_w(:, 2, :) = (pk_w(:, 2, :) - 1) / self.zero_padding_ratio;
                if length(intersect(round(pk_u(:,2,:)), round(pk_w(:,2,:)))) < self.num_user && ...
                        length(intersect(floor(pk_u(:,2,:)), floor(pk_w(:,2,:)))) < self.num_user
                    error('Error. Co-located preamble!'); 
                end
                for ui = 1: self.num_user
                    if std(pk_u(ui, 2, :)) >= 0.5
                        error('Error. Co-located preamble!'); 
                    end
                end
            end
            pk_u(pk_u<0) = [];
            % sort peak with index
            pk_u_i = pk_u(:, 2, :);
            pk_u_i(pk_u_i>0.5*self.num_bb_samples) = pk_u_i(pk_u_i>0.5*self.num_bb_samples) - self.num_bb_samples;
            pk_u(:, 2, :) = pk_u_i;
            for ii = 1: self.preamble_len
               pk_u(:, :, ii) = sortrows(pk_u(:, :, ii), 2); 
            end
            pk_u_h = mean(pk_u(:, 1, :), 3);
            pk_u_i = mean(pk_u(:, 2, :), 3);
            
            % preamble sfd shift
            pk_d = self.dechirp((self.preamble_len+2)*self.num_samples+1, 2, self.upchirp_cw);
            pk_d(:, 2, :) = (pk_d(:, 2, :) - 1) / self.zero_padding_ratio;
            if self.debug_ppa
                pk_w = self.dechirp((self.preamble_len+2)*self.num_samples+1, 2, self.upchirp_cw);
                pk_w(:, 2, :) = (pk_w(:, 2, :) - 1) / self.zero_padding_ratio;
                if length(intersect(round(pk_d(:,2,:)), round(pk_w(:,2,:)))) < self.num_user && ...
                        length(intersect(floor(pk_d(:,2,:)), floor(pk_w(:,2,:)))) < self.num_user
                    error('Error. Co-located preamble!'); 
                end
                for ui = 1: self.num_user
                    if std(pk_d(ui, 2, :)) >= 0.5
                        error('Error. Co-located preamble!'); 
                    end
                end
            end
            pk_d(pk_d<0) = [];
            % sort peak with index
            pk_d_i = pk_d(:, 2, :);
            pk_d_i(pk_d_i>0.5*self.num_bb_samples) = pk_d_i(pk_d_i>0.5*self.num_bb_samples) - self.num_bb_samples;
            pk_d(:, 2, :) = pk_d_i;
            pk_d(:, :, 1) = sortrows(pk_d(:, :, 1), 2);
            pk_d(:, :, 2) = sortrows(pk_d(:, :, 2), 2);
            pk_d_h = mean(pk_d(:, 1, :), 3);
            pk_d_i = mean(pk_d(:, 2, :), 3);
            
            % co-located preamble
            if self.debug_ppa
                if size(pk_u, 1) < self.num_user || size(pk_d, 1) < self.num_user
                   error('Error. Co-located preamble!'); 
                end
                if ~isempty(find(abs(sort(pk_u_h)-sort(pk_d_h)) >= 30, 1))
                    error('Error. Co-located preamble!'); 
                end
            end
            
            % coarse cfo & to estimation
            % match cfo&to with peak magnitude, factorial(n) = n!
            n_fac = factorial(self.num_user);
            cfo_to_u = zeros(self.num_user, 2, n_fac);
            match_err = zeros(n_fac, 1);
            % n! permutations to combine pk_d with pk_u
            pk_d_i_perm = pk_d_i(perms(1: self.num_user)).';
            pk_d_h_perm = pk_d_h(perms(1: self.num_user)).';
            for fi = 1: n_fac
               % cfo
               cfo_to_u(:, 1, fi) = (pk_u_i + pk_d_i_perm(:, fi)) / 2;
               % to
               cfo_to_u(:, 2, fi) = (pk_u_i - pk_d_i_perm(:, fi)) / 2;
               % mag
               match_err(fi) = sum(abs(pk_u_h - pk_d_h_perm(:, fi)));
            end
            [~, min_i] = min(match_err);
            cfo_coarse = cfo_to_u(:, 1, min_i);
            to_coarse = cfo_to_u(:, 2, min_i);
            % assume that the signal always arrive after the demodulation window 
            if ~isempty(find(to_coarse > 0, 1))
               warning("Error. Invalid demodulation window."); 
               to_coarse = -0.1*ones(self.num_user, 1);
            end
            

            sig_ft = fft(self.sig(1: 2: self.num_samples) .* self.downchirp(1: 2: end));
            sig_e_t = self.sig_constr("time", zeros(self.num_user, 1), cfo_coarse, to_coarse, true);
            sig_e_ft = fft(sig_e_t(1: 2: end, :) .* self.downchirp(1: 2: end));
            chan_coarse = LoRaPPA.channel_lse(sig_e_ft, sig_ft);
            if self.debug_ppa
               figure(101);plot(abs(sig_ft));hold on;plot(abs(sig_e_ft));hold off; 
            end
            self.cfo_coarse_u = cfo_coarse;
            self.to_coarse_u = to_coarse;
            self.chan_coarse_u = chan_coarse;
            
            if self.pream_est == "fine"
                % joint optimization of cfo & to & chan
                % FIXME: looks like that fine optimization has little gain
                fft_bin_size = 1/self.zero_padding_ratio;
                cfo_ = zeros(self.num_user, 1, self.preamble_len);
                to_ = zeros(self.num_user, 1, self.preamble_len);
                chan_ = zeros(self.num_user, 1, self.preamble_len);
                for ui = 1: self.preamble_len
                    lb = [cfo_coarse-fft_bin_size to_coarse+fft_bin_size -inf*ones(self.num_user,1) -inf*ones(self.num_user,1)];
                    ub = [cfo_coarse+fft_bin_size to_coarse-fft_bin_size inf*ones(self.num_user,1) inf*ones(self.num_user,1)];
                    [opt_sol, opt_fval] = self.residual_min("offset+channel", "constrained", ft_u(:, ui), "freq", zeros(self.num_user, 1), cfo_coarse, to_coarse, true, (ui-1)*self.num_bb_samples, chan_coarse, lb, ub, []);
                    cfo_(:,:,ui) = opt_sol(:, 1);
                    to_(:,:,ui) = opt_sol(:, 2);
                    chan_(:,:,ui) = opt_sol(:, 3) + 1i*opt_sol(:, 4);
                end
                self.cfo_u = mean(cfo_, 3);
                self.to_u = mean(to_, 3);
                self.chan_u = mean(chan_, 3);
            elseif self.pream_est == "coarse"
                self.cfo_u = cfo_coarse;
                self.to_u = to_coarse;
                self.chan_u = chan_coarse;
            else
                error('Error. Invalid pream_est.');
            end
            
            % sfo estimation
            if self.is_real
                self.sfo_drift = ((1 + (1:self.symbol_len)') * 2^self.sf * self.cfo_u.'*self.bw / (self.rf_freq*self.num_bb_samples)).'; % row: user_idx, col: syn_idx
            else
                self.sfo_drift = zeros(self.num_user, self.symbol_len);
            end
            
            if self.debug_ppa
                sig_e = self.sig_constr("freq", zeros(self.num_user), self.cfo_u, self.to_u, true, 0, self.chan_u);
                figure(101);
                plot(ft_u(:, 1), '.-');
                hold on;
                plot(sum(sig_e, 2), '.-');
                hold off;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             global amp;
%             global phase;
%             global cfo;
%             global to; 
%             self.chan_u = amp .* phase;
%             self.cfo_u = cfo*self.num_bb_samples/self.bw;
%             self.to_u = -to;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
        function y = sig_constr(self, domain, sym, cfo, to, is_up, cfo_start, chan)
            % sig_constr  reconstruct chirps
            %
            % input:
            %     domain: reconstruction domain, `time` or `freq`
            %     sym: a column vector of symbols
            %     cfo: cfo (BIN) vector of users
            %     to: to vector of users
            %     is_up: `true` if applying up-chirp dechirping
            %            `false` if applying down-chirp dechirping
            %     cfo_start: start index of cfo
            %     chan: signal channel
            % 
            % output:
            %     y: reconstructed chirp
            %
            
            if nargin == 6
                cfo_start = 0;
                chan = ones(size(sym));
            elseif nargin == 7
                chan = ones(size(sym));
            end

            sym = sym(:, 1);
            sym = mod(sym, self.num_bb_samples);
            num_sym = length(sym);
            fs_bw = self.fs_constr / self.bw;
            if is_up
                c = self.downchirp_c;
            else
                c = self.upchirp_c;
            end
            y_ = zeros(self.num_samples, num_sym);
            to = round(abs(fs_bw*to));
            for si = 1: num_sym
                if self.fast_mode && self.algo ~= "unc"
                    sig = self.fast_chirp(:, sym(si) + 1);
                else
                    sig = LoRaPPA.chirp(self.sf, self.bw, self.fs_constr, sym(si), is_up);
                end
                % add cfo shift
                sig = LoRaPPA.add_cfo(sig, cfo(si), self.num_bb_samples, self.num_constr_samples, cfo_start);
                % add to
                sig = [zeros(to(si), 1); sig];
                sig = chan(si) * sig(1: fs_bw*self.num_bb_samples);
                y_(:, si) = sig(1: fs_bw/2: end);
            end
            
            if domain == "freq"
                y = abs(fft(sum(y_,2) .* c, self.fft_size));
                y = y(1:self.bin_size, :) + y(self.fft_size-self.bin_size+1:self.fft_size, :);
            elseif domain == "time"
                y = y_;
            end
        end
        
        function [data, crc_check] = soft_decode(self, symbols_p0)
            % soft_decode  Decode symbols to bytes
            %
            % input:
            %     symbols_pr0: A matrix containing probabilities of
            %                demodulated bits that equal to zero.
            % output:
            %     data: Decoded payload of LoRa packet
            %     crc_check: valid: 1, invalid: 0
            
            % soft gray encode
            symbols_g_bin = self.soft_gray_encoding(symbols_p0);
            % bit-level deinterleave
            symbols_di_bin = self.deinterleave_bin(symbols_g_bin);
            % soft hamming decode
            nibbles = self.soft_hamming_decode(symbols_di_bin);
            % nibbles to bytes
            bytes = uint8(zeros(min(255, floor(length(nibbles)/2)), 1));
            for ii = 1:length(bytes)
                bytes(ii) = bitor(uint8(nibbles(2*ii-1)), 16*uint8(nibbles(2*ii)));
            end
            % dewhitening
            plen = self.payload_len;
            if self.crc
                % last 2 bytes are CRC16 checkcum
                data = self.dewhitening(bytes(1:plen));
                % calculate CRC checksum
                checksum = self.calc_crc(data(1:plen));
                if isequal(checksum, bytes(plen+1:plen+2))
                    crc_check = true;
                else
                    crc_check = false;
                end
            else
                data = self.dewhiten(bytes(1:len));
                crc_check = false;
            end
        end
        
        function symbols_g_bin = soft_gray_encoding(self, symbols_p0)
            % soft_gray_encoding  Soft gray encoding 
            %
            % input:
            %     symbols_pr0: A matrix containing probabilities of
            %                demodulated bits that equal to zero.
            % output:
            %     symbols_g_bin: soft gray encoding bits
            % 
            
            % vector to matrix
            symbols_p0_ = reshape(symbols_p0, [self.sf, self.symbol_len]).';
            % do the shift of -1
            symbols_p0_m1 = zeros(self.symbol_len, self.sf);
            for si = 1: self.symbol_len
                symbols_p0_m1(si, :) = LoRaPPA.soft_bits_m1(symbols_p0_(si, :));
            end
            % B >> 1
            symbols_p0_sr1 = symbols_p0_m1(:, 1: self.sf-1);
            symbols_p0_sr1 = [ones(self.symbol_len, 1) symbols_p0_sr1];
            % C = B XOR B>>1
            symbols_g_bin = zeros(self.symbol_len, self.sf);
            for si = 1: self.symbol_len
                for bi = 1: self.sf
                    symbols_g_bin(si, bi) = LoRaPPA.soft_bits_xor(symbols_p0_m1(si, bi), symbols_p0_sr1(si, bi));
                end
            end
        end
        
        function nibbles = soft_hamming_decode(self, codewords_p0)
            % soft_hamming_decode  Soft sphere hamming decoding
            % ref to https://jp.mathworks.com/matlabcentral/fileexchange/42953-soft-hamming-decoder?s_tid=FX_rc2_behav
            % FIXME: how to convert probability to LLR
            % nsdec = 5;      
            % symbols_ll = (codewords_p0 - 0.5) * 2^nsdec;
            
            nsdec = 4;
            symbols_ll = log(codewords_p0./(1-codewords_p0));
            symbols_ll(symbols_ll>2^nsdec) = 2^nsdec;
            symbols_ll(symbols_ll<-2^nsdec) = -2^nsdec;
            
            nibbles_m = [];
            cw_len = 4;
            % use reduced rate for the first block
            cr_first = 4;
            sf_app = self.sf - 2;
            [HT, e_sort, s_pdf] = LoRaPPA.gen_hamming(cr_first);
            [x_bit, x_ll] = LoRaPPA.hamm_soft_out(HT, e_sort, s_pdf, symbols_ll(1:sf_app, :));
            nibbles_m = [nibbles_m; fliplr(x_bit(:, 1:cw_len))];
            
            % Normal rate for following block
            [HT, e_sort, s_pdf] = LoRaPPA.gen_hamming(self.cr);
            [x_bit, x_ll] = LoRaPPA.hamm_soft_out(HT, e_sort, s_pdf, symbols_ll(sf_app+1:end, 4-self.cr+1:end));
            nibbles_m = [nibbles_m; fliplr(x_bit(:, 1:cw_len))];
            nibbles = [];
            for ni = 1: size(nibbles_m, 1)
                nibbles = [nibbles bi2de(nibbles_m(ni, :), 'left-msb')];
            end
        end
        
        function y = peak_enum(self, pk_i, m, algo, pk_h, th)
            % peak_enum  peak enumeration with given row vector and number
            %
            % input:
            %     pk_i: peak index vector
            %     m: number of users
            %     algo: enumeration algorithm
            %     pk_h: peak height vector
            %     th: peak match threshold for `Match-M!` algo
            % 
            % output:
            %     y: generated peak permutations
            %
            if nargin == 3
                algo = "M^M";
                pk_h = [];
                th = 0.1;
            elseif nargin == 4
                pk_h = [];
                th = 0.1;
            elseif nargin == 5
                th = 0.1;
            end
            if ~isrow(pk_i)
                pk_i = pk_i.';
            end
            if algo == "Match-M!"
                chan_h = abs(self.chan_u);
                assert(length(pk_i) == m && length(pk_i) == length(pk_h) && issortedrows(pk_h, 'descend') && issortedrows(chan_h, 'descend'));
                % peak match first
                chan_h = abs(self.chan_u);
                pk_h = pk_h / self.num_bb_samples;
                pk_match = zeros(m, 1);     % indicate if the peak is matched
                num_match = 0;
                for nn = 1: m
                    hi = find(abs(pk_h - chan_h(nn)) < th, 1);
                    if hi
                       pk_match(hi) = 1; 
                       num_match = num_match + 1;
                    end
                end
                % enumeration
                if num_match
                    % all peaks match
                    if num_match == m
                        y = pk_i;
                        return;
                    end
                    % part of peaks match
                    y_ = self.peak_enum(pk_i(pk_match < 1), m - num_match, "M!");
                    len_perm = size(y_, 1);
                    y = [];
                    for nn = 1: m
                        if pk_match(nn)
                            y = [y pk_i(nn) * ones(len_perm, 1)];
                        else
                            y = [y y_(:, 1)];
                            y_(:, 1) = [];
                        end
                    end
                else
                    % no peak matches
                    y = self.peak_enum(pk_i, m, "M!");
                end
            elseif algo == "M!"
                assert(length(pk_i) == m);
                y = [];
                for nn = 1: length(pk_i)
                    y = [y; LoRaPPA.vec_perm(pk_i(1:nn), m)];
                end
            elseif algo == "M^M"
                assert(length(pk_i) == m);
                % permutation with repetition
                vec = LoRaPPA.nmultichoosek(pk_i, m);
                y = [];
                for i = 1: size(vec, 1)
                    y = [y; perms(vec(i,:))];
                end
                y = unique(y, 'rows');
            elseif algo == "V^M"
                assert(length(pk_i) >= m);
                if length(pk_i) > m
                    y = LoRaPPA.vec_perm(pk_i, m);
                else
                    y = self.peak_enum(pk_i, m, "M^M");
                end
            else
               error("Invalid algo!") 
            end
        end
        
        function [pk, ft_] = dechirp(self, x, n, chirp)
            % dechirp  Apply dechirping on the multiple symbols
            %
            % input:
            %     x: Start index of symbols
            %     n: number of symbols to be dechirped
            %     chirp: chirp to be applied
            % 
            % output:
            %     pk: Peak in FFT results of dechirping
            %         pk = [(height_1, index_1);
            %               (height_2, index_2)]
            %     ft_: fft magnitude summation
            %

            c = repmat(chirp, n, 1);
            ft = fft(reshape(self.sig(x:x+n*self.num_samples-1).*c, self.num_samples, n), self.fft_size);
            ft_ = abs(ft(1:self.bin_size, :)) + abs(ft(self.fft_size-self.bin_size+1:self.fft_size, :));
            pk = zeros(self.num_user, 2, n);
            for nn = 1: n
                pk_ = self.top_pks(ft_(:, nn), self.num_user);
                % padding -1 to match size
                pk(:, :, nn) = [pk_; -ones(self.num_user - size(pk_,1), 2)];
                %self.plot_fft(ft_(:, nn), pk_);
            end
        end
        
        function y = top_pks(self, pks, n, th)
            if nargin == 3
                th = self.threshold;
            end
            pks = abs(pks);
            [~, l] = findpeaks(pks);
            
            % window exchange to find peaks at boundary
            mid = round(self.bin_size/2);
            pks_exc = [pks(mid+1: end); pks(1: mid)];
            [~, l_exc] = findpeaks(pks_exc);
            l_exc_l = l_exc(l_exc > mid) - mid;
            l_exc_r = l_exc(l_exc <= mid) + mid;
            % pks union
            l = union(l, [l_exc_l; l_exc_r]);
            val = pks(l);
            % higher than th
            idx = find(val > th);
            l = l(idx);
            val = val(idx);
            % topn
            [~, n_idx] = maxk(val, n);
            l = l(n_idx);
            val = pks(l);
            y = [val l];
        end
    end
    
    methods(Static)
        function [corr_norm, corr_val] = cross_correlate(raw_samples, start_index, known)
            % cross_correlate  time-domain cross correlation
            % A = symbol[1:M]
            % B = d_known[1:M]
            % corr = abs(A .* conj(B)) / sqrt(norm(A)*norm(B))
            % 
            % input:
            %     raw_samples: raw samples from hardware
            %     start_index: correlation start index
            %     known: known time-domain signals to correlate with
            % 
            % output:
            %     corr_norm: normalized correlation values
            %     corr_val: correlation values
            %  
            p_known = sum(known .* conj(known));
            window_len = length(known);
            corr_norm = [];
            corr_val = [];
            for i = start_index: start_index + length(raw_samples) - 1 - window_len  % may user PARFOR here
                corr = sum(raw_samples(i: i + window_len - 1) .* conj(known));
                p = sum(raw_samples(i: i + window_len - 1) .* conj(raw_samples(i: i + window_len - 1)));
                corr_norm = [corr_norm; abs(corr) / sqrt(p * p_known)];
                corr_val = [corr_val; abs(corr)];
            end
        end 
       
        function h = channel_lse(e, sig)
            % channel_lse  least-square estimation of channel
            %
            % input:
            %     e: estimation signal
            %     sig: rx signal
            % 
            % output:
            %     h: lse of channel
            %
            
           h = inv((e.') * e) * (e.') * sig;
        end
        

        
        function r = residual_func_bin(sig, sig_e, domain, bins)
            % residual_func_bin  residual function of top peak bins
            %
            % input:
            %     sig: rx signal
            %     sig_e: estimation signal
            %     domain: freq domain reconstruction signal
            % 
            % output:
            %     r: residual value
            %
            assert(domain == "freq");
            r = norm(abs(sig(bins)) - abs(sig_e(bins)), 2)^2;
        end
        
        function y = vec_perm(n, m)
            % vec_perm  generate permuatitions with given vector and number
            len_n = length(n);
            ex_len = m - len_n;       
            if ex_len > 0                              % M>N 
                vec_extra = LoRaPPA.nmultichoosek(n, ex_len);   % repeatedly extra in N
                len_vec_ex = size(vec_extra, 1);
                vec = zeros(len_vec_ex, m);
                vec(:, 1:len_n) = repmat(n, len_vec_ex, 1);  % former N is identical
                vec(:, len_n+1:end) = vec_extra;           % all possible combinations
                y=[];
                for i = 1: len_vec_ex
                    y = [y; perms(vec(i,:))];              % all possible permutations
                end
                y = unique(y, 'rows');                     % remove duplicate
            elseif ex_len < 0                         % M<N
                % permutation with repetition
                vec = LoRaPPA.nmultichoosek(n, m);
                y = [];
                for i = 1: size(vec,1)
                    y = [y; perms(vec(i,:))];
                end
                y = unique(y, 'rows');
            else
                y = perms(n);
            end   
        end

        function combs = nmultichoosek(values, k)
            if numel(values)==1 
                n = values;
                combs = nchoosek(n+k-1,k);
            else
                n = numel(values);
                combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
                combs = reshape(values(combs),[],k);
            end
        end
        
        function y = soft_bits_m1(p0, flg)
            % soft_bits_m1  do shift of -1 to bits with p0
            % prev_pr_0 = prod(pr_0(curr-1:end))
            % out_pr_0 = curr_pr_1*prev_pr_0 + curr_pr_0*(1-prev_pr_0)
            % e.g. soft_bits_minus_1([0.4, 0.2, 0.7])->[0.428, 0.62, 0.3]
            if nargin == 1
                flg = "left-msb";
            end
            if flg == "left-msb"
                len = length(p0);
                y = zeros(1, len);
                % flip for last bit
                y(len) = 1 - p0(len);
                % transfer for other bits
                for bi = 1: len -1
                    prev_p0 = prod(p0(bi+1: len));
                    y(bi) = (1 - p0(bi)) * prev_p0 + p0(bi) * (1 - prev_p0);
                end
            elseif flg == "right-msb"
                error('Error. Not support right-msb now.');
            else
                error('Error. Invalid flag.');
            end
        end
        
        function y = soft_bits_xor(b1_p0, b2_p0)
            % soft_bits_xor  XOR of two soft bits with p0
            y = b1_p0 * b2_p0 + (1 - b1_p0) * (1 - b2_p0);
        end
        
        function [x_decode,x_llr_new, count] = hamm_soft_out(HT, e_sort, s_pdf, x_llr)
            % hamm_soft_out  calculates the soft-output x_llr_new from the soft-input x_llr

            Nframes=size(x_llr,1);
            nc=size(x_llr,2);

            x_llr(x_llr==0)=0.000000001;

            % Memory allocation
            x_p_new=NaN(1,nc);
            x_llr_new=x_llr;

            x_p0=exp(abs(x_llr))./(1+exp(abs(x_llr)));

            x_code_r=(sign(-x_llr)+1)/2;
            x_decode=x_code_r;
            count=0;
            % Decoding, old version (tested, slow)
            x_decode=x_code_r;
            for k=1:Nframes
                % Syndrom calculation
                syn_b=mod(x_code_r(k,:)*HT,2);
                if sum(syn_b)~=0
                    count=1;
                    %syn_d=bi2de(syn_b,[],'left-msb'); % too slow
                    syn_b = syn_b(end:-1:1);
                    pow2vector = 2.^(0:1:(size(syn_b,2)-1));
                    syn_d = syn_b*pow2vector';

                    %find start and end in s_sort for syn_b
                    start=sum(s_pdf(1:syn_d-1))+1;
                    ende=sum(s_pdf(1:syn_d));
                    stepsize=ende-start+1; % number of error patterns for syndrome
                    e_work=e_sort(start:ende,:);

                    %multiply row-wise and find minimum, decode to minimum
                    e_llr_prod=LoRaPPA.bem_repmat(abs(x_llr(k,:)),stepsize).*e_work;
                    [~,idx_llr_min]=min(sum(abs(e_llr_prod),2));
                    x_decode(k,:)=mod(x_code_r(k,:)+e_sort(start+idx_llr_min-1,:),2);

                    % x_p0 probabilit? que bit soit juste
                    % x_p probabilit? que bit (flipp?) soit juste, pour chaque pattern
                    x_p=abs(e_work-LoRaPPA.bem_repmat(x_p0(k,:),stepsize));
                    % P: probabilit? qu'un des errors patterns ait eu lieu
                    P=prod(x_p,2).';

                    % Normalize Probability Pnorm
                    Psum=sum(P);
                    Pnorm=P./Psum;
                    if Psum==0
                        Pnorm(:)=0.000001;
                    end

                    % Limit too high propabilities
                    if max(Psum)>(1-1E-25)
                        Pnorm=Pnorm*(1-1E-25);
                    end

                    for k4=1:nc % bit par bit du codeword
                        wert=mod(x_code_r(k,k4)+e_work(:,k4),2);
                        x_p_new(k4)=1-sum(wert.*Pnorm.');
                    end
                    % calculate new LLRs
                    x_llr_new(k,:)=log(x_p_new./(1-x_p_new));
                end


                % Limit infinit LLRs max of input
                Threshold=1*max(max(abs(x_llr)))+randn(1)/10;
                x_llr_new(k,x_llr_new(k,:)>Threshold)=Threshold;
                Threshold=Threshold+randn(1)/10;
                x_llr_new(k,x_llr_new(k,:)<-Threshold)=-Threshold;
            end
            x_llr_new=real(x_llr_new);
        end
        
        function [HT, e_sort, s_pdf] = gen_hamming(CR)
            % gen_hamm_params  Generation of Hamming matrices G and H
            if CR == 1
                G = [1,0,0,0,1
                     0,1,0,0,1
                     0,0,1,0,1
                     0,0,0,1,1];
                H = [1,1,1,1,1];
                HT = H.';
            else
                G = [1,0,0,0,1,0,1,1
                     0,1,0,0,1,1,1,0
                     0,0,1,0,1,1,0,1
                     0,0,0,1,0,1,1,1];
                G = G(:, 1:CR+4);
                H = gen2par(G);             % generate H from G
                HT = H.';                   % Transpose
            end
            % Preparing error patterns for soft-decoding
            % Permutation of 1,2 and 3 bit errors
            % FIXME: seems that 2 is the best
            [e1, e2, e3] = LoRaPPA.permute_e(4 + CR, 2);
            e = [e1;e2;e3];

            % syndrom calculation
            syn = NaN(size(e,1), CR);
            for k = 1:size(e,1)
                syn(k,:) = mod(e(k,:)*HT, 2);
            end

            % sort rows of syndroms associated error pattern e
            [s_sort,s_idx] = sortrows(syn);     % sort syn to s_sort
            e_sort=e(s_idx,:);                  % sort e to e_sort in the same way as syn

            e_sort=e_sort(sum(s_sort,2)>0,:);   % delete zero rows
            s_sort=s_sort(sum(s_sort,2)>0,:);   % delete zero rows

            s_de=LoRaPPA.bem_bi2de(s_sort);             % calculate decimal values from s_sort

            s_pdf=hist(s_de,1:2^CR-1);          % calculate which value occurs how often
        end
        
        function [eM1, eM2, eM3] = permute_e(nc, M)
            %generates the bit error pattern e for a given codewordlength nc
            %up to the maximium error number M

            % Permutation of 1 bit error
            eM1 = logical(eye(nc));

            % Permutation of 2 bit error
            eM2 = [];
            if M > 1
                %MeM2=round(factorial(nc)/(factorial(nc-2)*factorial(2)));
                MeM2 = (nc*(nc-1))/2;
                eM2 = [1 1 zeros(1,nc-2); NaN(MeM2-1,nc)];
                z = 1;
                for ko = 0: nc - 3
                    for ki = 2: nc - 1 - ko
                        z = z + 1;
                        eM2(z, :) = 0;
                        eM2(z, ko+1) = 1;
                        eM2(z, ki+ko)=0;
                        eM2(z, ki+ko+1)=1;
                    end
                    z = z + 1;
                    eM2(z, :) = 0;
                    eM2(z, ko+2) = 1;
                    eM2(z, ko+3) = 1;
                end
                eM2 = logical(eM2);
            end

            % Permutation of 3 bit error
            eM3 = [];
            if M > 2
                % MeM3=factorial(nc)/(factorial(nc-3)*factorial(3));
                MeM3 = (nc*(nc-1)*(nc-2)) / 6;
                eM3 = false(MeM3, nc);
                arg = 0;
                for ko = 0: nc-3
                    for km = 0: nc-3
                        if nc-2-km-ko<1; break; end
                        e_sub = logical([zeros(nc-2-km-ko,ko) ones(nc-2-km-ko,1)...
                            zeros(nc-2-km-ko,km) ones(nc-2-km-ko,1), eye(nc-2-km-ko)]);
                        arg = (1:size(e_sub,1)) + max(arg);
                        eM3(arg,:) = logical(e_sub);
                    end
                end
            end
        end
        
        function d=bem_bi2de(b)
            % simple function for calculating decimal d form binary b

            b = b(:,end:-1:1); %'left-msb'
            pow2vector = 2.^(0:1:(size(b,2)-1));
            size_B = size(b,2);
            d = b(:,1:size_B)*pow2vector(:,1:size_B).';
        end
        
        function B = bem_repmat(A,M)
            % simple repmat function 
            % generates a matrix B with M row out of a column vector A

            B = A(ones(M, 1), :);
        end
    end
end