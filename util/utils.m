classdef utils < handle
    properties
        color_blue
        color_orange
        color_yellow
        color_purple
        color_green
        color_cyan
        color_red
        color_black
    end
    
    methods
        function self = utils()
            self.color_blue = [0 0.4470 0.7410];
            self.color_orange = [0.8500 0.3250 0.0980];
            self.color_yellow = [0.9290 0.6940 0.1250];
            self.color_purple = [0.4940 0.1840 0.5560];
            self.color_green = [0.4660 0.6740 0.1880];
            self.color_cyan = [0.3010 0.7450 0.9330];
            self.color_red = [0.6350 0.0780 0.1840];
            self.color_black = [0 0 0];
        end
    end
        
    methods (Static)

        function dataout = add_noise_trace_driven(datain, snr, data_noise, amp_sig)
            % snr = 10*log((s^2/n_0^2)/n^2)
            if nargin == 3
                amp_sig = mean(abs(datain));
            end
            amp_noise_level = mean(abs(data_noise));
            amp_noise = amp_sig / 10^(snr/20) / amp_noise_level;
            dataout  = datain + amp_noise/sqrt(2) * data_noise;
        end
        
        function dataout = add_noise(datain, snr)
            % snr = 10*log(s^2/n^2)
            amp_sig = mean(abs(datain));
            amp_noise = amp_sig / 10^(snr/20);
            dlen = length(datain);
            dataout  = datain + (amp_noise/sqrt(2) * randn([dlen 1]) + 1i*amp_noise/sqrt(2) * randn([dlen 1]));
        end
        
        function plot_curve(fn, x, y, options)
            arguments
               fn;
               x;
               y;
               options.Marker = ['o','s','d','^','v','x','*','+'];
               % blue, orange, yellow, purple, green, cyan, red, black
               options.Color = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560];...
                                [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]; [0 0 0]];
               options.MarkerSize = 5;
               options.LineWidth = 1.8;
            end
            figure;
            % fn: plot/semilogx/semilogy
            marker = options.Marker;
            color = options.Color;
            assert(length(x) == size(y, 2));
            for i = 1: size(y, 1)
               fn(x, y(i, :), [marker(i) '-'], 'Color', color(i, :), 'MarkerSize', options.MarkerSize, 'LineWidth', options.LineWidth);
               hold on;
            end
            hold off;
            set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5);
            grid on;
        end
        
        function [res, idx] = perm_match(data, base)
           % match base
           perm_data = perms(data).';
           perm_idx = perms(1: length(data)).';
           n_fac = factorial(length(data));
           err = zeros(n_fac, 1);
           for ii = 1: n_fac
                err(ii) = sum(abs(perm_data(:,ii) -  base), 'all');
           end
           [~, min_i] = min(err);
           idx = perm_idx(:, min_i);
           res = data(idx);
        end
        
        function [rat, num] = calc_ser(data, gt)
            % sort for calculation
            len = size(gt, 1);
            n = size(gt, 2);
            n_fac = factorial(n);
            perm_ = perms(1: n);
            err = zeros(n_fac, 1);
            for ii = 1: n_fac
                perm_gt = gt(:, perm_(ii,:));
                err(ii) = sum(abs(data - perm_gt), 'all');
            end
            [~, min_i] = min(err);
            gt = gt(:, perm_(min_i,:));

            num = 0;
            for ui = 1: n
                [~, num_] = utils.calc_ser_(data(:,ui), gt(:,ui));
                num = num + num_;
            end
            rat = num/(len*n);
        end

        function [rat, num] = calc_ser_(data, gt)
            len1 = length(data);
            len2 = length(gt);
            len = min(len1, len2);
            if len > 0
                data = data(1: len);
                gt = gt(1: len);
                num = sum(data(:) ~= gt(:), 'all') + (len2 - len);
                rat = num / len2;
            else
                num = len2;
                rat = 1;
            end
        end
        
        function [rat, num, perm_gt] = calc_ber(data, gt)
            % sort for calculation
            len = size(gt, 1);
            n = size(gt, 2);
            n_fac = factorial(n);
            perm_ = perms(1: n);
            err = zeros(n_fac, 1);
            for ii = 1: n_fac
                perm_gt = gt(:, perm_(ii,:));
                err(ii) = sum(abs(double(data) - double(perm_gt)), 'all');
            end
            [~, min_i] = min(err);
            perm_gt = perm_(min_i,:);
            gt = gt(:, perm_gt);

            num = 0;
            for ui = 1: n
                [~, num_] = utils.calc_ber_(data(:,ui), gt(:,ui));
                num = num + num_;
            end
            rat = num/(len*n*8);
        end

        function [rat, num] = calc_ber_(data, gt)
            len1 = length(data);
            len2 = length(gt);
            len = min(len1, len2);
            if len > 0
                data = data(1: len);
                gt = gt(1: len);
                num = sum(xor(de2bi(data, 8, 2, 'left-msb'), de2bi(gt, 8, 2, 'left-msb')), 'all') + 8 * (len2 - len);
                rat = num / (len2 * 8);
            else
                num = len2 * 8;
                rat = 0.5;
            end
        end
        
        function [rat, num] = calc_sum_ber(data, gt)
            assert(isequal(size(data), size(gt)));
            [len, n] = size(gt);
            data_sum = sum(reshape(de2bi(data, 8, 2, 'left-msb').', 8*len, n), 2);
            gt_sum = sum(reshape(de2bi(gt, 8, 2, 'left-msb').', 8*len, n), 2);
            num = sum(data_sum ~= gt_sum);
            rat = num / (8*len);
        end
        
        % Method `read` and `write` are copied from
        % https://github.com/gnuradio/gnuradio/blob/master/gr-utils/octave/read_complex_binary.m
        % https://github.com/gnuradio/gnuradio/blob/master/gr-utils/octave/write_complex_binary.m
        function v = read(filename, count)

            % usage: read(filename, [count])
            %
            %  open filename and return the contents as a column vector,
            %  treating them as 32 bit complex numbers
            %

            m = nargchk (1,2,nargin);
            if (m)
            usage (m);
            end

            if (nargin < 2)
            count = Inf;
            end

            f = fopen (filename, 'rb');
            if (f < 0)
            v = 0;
            else
            t = fread (f, [2, count], 'float');
            fclose (f);
            v = t(1,:) + t(2,:)*1i;
            [r, c] = size (v);
            v = reshape (v, c, r);
            end
        end

        function v = write(data, filename)

            % usage: write(data, filename)
            %
            %  open filename and write data to it
            %  Format is interleaved float IQ e.g. each
            %  I,Q 32-bit float IQIQIQ....
            %  This is compatible with read_complex_binary()
            %

            m = nargchk (2,2,nargin);
            if (m)
            usage (m);
            end

            f = fopen (filename, 'wb');
            if (f < 0)
            v = 0;
            else
            re = real(data);
            im = imag(data);
            re = re(:)';
            im = im(:)';
            y = [re;im];
            y = y(:);
            v = fwrite (f, y, 'float');
            fclose (f);
            end
        end
    end
end