clear; clc
% Parameter
FFT_Size = 128;
FFT_Size_SSKM = FFT_Size * 2;

GI_Size = FFT_Size / 4;
GI_Size_SSKM = FFT_Size_SSKM / 4;

Modulation_Order = 4;

Data_Size = FFT_Size * Modulation_Order;
Data_Size_SSKM = Data_Size * 2;

Multi_Path = 7;
SNR = 0:3:30;
Iteration = 10000;

BER_SC = zeros(1, length(SNR));
BER_MRC = zeros(1, length(SNR));
BER1_SSKM = zeros(1, length(SNR));
BER2_SSKM = zeros(1, length(SNR));

bit_error_total_SC = 0;
bit_error_total_MRC = 0;
bit_error_total1_SSKM = 0;
bit_error_total2_SSKM = 0;

for i = 1:length(SNR)
    
    for j = 1:Iteration
        % 송신단 데이터 생성
        Data = randi([0 1], 1, Data_Size);

        
        % 변조 
        Modulate = base_mod(Data, Modulation_Order);

        % IFFT
        mod = ifft(Modulate, [], 2) * sqrt(FFT_Size);

        % Insert cp
        cp = mod(FFT_Size - GI_Size + 1:end);
        x = cat(2, cp, mod);

        % 채널 생성(channel, noise 고려)
        h1 = Rayleigh_channel(Multi_Path);
        h2 = Rayleigh_channel(Multi_Path);

        % 컨벌루션
        hx1 = conv(h1, x);
        hx2 = conv(h2, x);

        % 노이즈 첨가
        y1 = awgn(hx1, SNR(i));
        y2 = awgn(hx2, SNR(i));

        % CP 제거
        y_cp_removed1 = y1(GI_Size + 1:GI_Size + FFT_Size); % 33~160 까지 총 128개
        y_cp_removed2 = y2(GI_Size + 1:GI_Size + FFT_Size);

        % FFT
        Y_fft1 = fft(y_cp_removed1) / sqrt(FFT_Size)/sqrt(2);
        Y_fft2 = fft(y_cp_removed2) / sqrt(FFT_Size)/sqrt(2);

        H1 = fft(h1, FFT_Size);
        H2 = fft(h2, FFT_Size);

        % SC 알고리즘(subcarrier마다 큰 h를 select해서 최적의 h를 select)
        H_temp = [H1; H2];
        Y_sc = [Y_fft1; Y_fft2];

        H_sc = zeros(1, FFT_Size);
        Y_selected = zeros(1, FFT_Size);

        for k = 1:FFT_Size
            [~, index] = max(abs(H_temp(:, k)).^2);
            H_sc(k) = H_temp(index, k);
            Y_selected(k) = Y_sc(index, k);
        end


        % MRC 알고리즘(subcarrier마다 각h 별로 가중치를 곱해서 다 더해서 모든 채널에 대해 크기가 고려되어서 diversity가 결정됨.)
        W1 = conj(H1);
        W2 = conj(H2);
        Y_MRC = (Y_fft1 .* W1 + Y_fft2 .* W2);
        H_mrc = abs(H1).^2 + abs(H2).^2;

        % Equalizer
        Y_equalized_sc = Y_selected ./ H_sc;
        Y_equalized_mrc = Y_MRC ./ H_mrc;

        % Demod
        DeMod_SC = base_demod(Y_equalized_sc, Modulation_Order);
        DeMod_MRC = base_demod(Y_equalized_mrc, Modulation_Order);

        % BER_total
        bit_error_numb1 = find(DeMod_SC ~= Data);
        bit_error_numb2 = find(DeMod_MRC ~= Data);
     
        bit_error_total_SC  = bit_error_total_SC + length(bit_error_numb1);
        bit_error_total_MRC = bit_error_total_MRC + length(bit_error_numb2);

        % SSKM Encoding
        % 16QAM
        % Data_SSKM = zeros(1, Data_Size_SSKM);
        % for n = 1:Data_Size/2
        %     Data_SSKM(4*n-3) = Data(2*n-1);
        %     Data_SSKM(4*n-2) = Data(2*n-1);
        %     Data_SSKM(4*n-1) = Data(2*n);
        %     Data_SSKM(4*n) = Data(2*n);
        % end
         %QPSK
        Data_SSKM = zeros(1,Data_Size_SSKM);
        for n=1 : Data_Size
            if Data(n)==1
                Data_SSKM(2*n-1) = 1;
                Data_SSKM(2*n) = 0;
        
            elseif Data(n)== 0
                Data_SSKM(2*n-1) = 0;
                Data_SSKM(2*n) = 1;
            end
        end
        % 변조 
        Modulate_SSKM = base_mod(Data_SSKM, Modulation_Order);

        % IFFT
        mod_SSKM = ifft(Modulate_SSKM, [], 2) * sqrt(FFT_Size_SSKM);

        % Insert cp
        cp_SSKM = mod_SSKM(FFT_Size_SSKM - GI_Size_SSKM + 1:end);
        x_SSKM = cat(2, cp_SSKM, mod_SSKM);

        
        % 컨벌루션
        hx1_SSKM = conv(h1, x_SSKM);
        hx2_SSKM = conv(h2, x_SSKM);

        % 노이즈 첨가
        y1_SSKM = awgn(hx1_SSKM, SNR(i));
        y2_SSKM = awgn(hx2_SSKM, SNR(i));

        % CP 제거
        y_cp_removed1_SSKM = y1_SSKM(GI_Size_SSKM + 1:GI_Size_SSKM + FFT_Size_SSKM);
        y_cp_removed2_SSKM = y2_SSKM(GI_Size_SSKM + 1:GI_Size_SSKM + FFT_Size_SSKM);

        % FFT
        Y_fft1_SSKM = fft(y_cp_removed1_SSKM) / sqrt(FFT_Size_SSKM);
        Y_fft2_SSKM = fft(y_cp_removed2_SSKM) / sqrt(FFT_Size_SSKM);

        H1_SSKM = fft(h1, FFT_Size_SSKM);
        H2_SSKM = fft(h2, FFT_Size_SSKM);

        % SC 알고리즘(subcarrier마다 큰 h를 select해서 최적의 h를 select)
        H_temp_SSKM = [H1_SSKM; H2_SSKM];
        Y_sc_SSKM = [Y_fft1_SSKM; Y_fft2_SSKM];

        H_sc_SSKM = zeros(1, FFT_Size_SSKM);
        Y_selected_SSKM = zeros(1, FFT_Size_SSKM);

        for k = 1:FFT_Size_SSKM
            [~, index] = max(abs(H_temp_SSKM(:, k)).^2);
            H_sc_SSKM(k) = H_temp_SSKM(index, k);
            Y_selected_SSKM(k) = Y_sc_SSKM(index, k);
        end

        % MRC 알고리즘
        W1_SSKM = conj(H1_SSKM);
        W2_SSKM = conj(H2_SSKM);
        Y_MRC_SSKM = (Y_fft1_SSKM .* W1_SSKM + Y_fft2_SSKM .* W2_SSKM)*sqrt(2);
        H_mrc_SSKM = abs(H1_SSKM).^2 + abs(H2_SSKM).^2; 

        % Equalizer
        Y_equalized_sc_SSKM = Y_selected_SSKM ./ H_sc_SSKM;
        Y_equalized_mrc_SSKM = Y_MRC_SSKM ./ H_mrc_SSKM;

        % Demod
        DeMod1_SSKM = base_demod(Y_equalized_sc_SSKM, Modulation_Order);
        DeMod2_SSKM = base_demod(Y_equalized_mrc_SSKM, Modulation_Order);

        %Dessk
        dessk_SC = deSSK(DeMod1_SSKM);
        dessk_MRC = deSSK(DeMod2_SSKM);

        %Bit error_total
        bit_error_numb1_SSKM = find(dessk_SC ~= Data);
        bit_error_numb2_SSKM = find(dessk_MRC ~= Data);
     
        bit_error_total1_SSKM = bit_error_total1_SSKM + length(bit_error_numb1_SSKM);
        bit_error_total2_SSKM = bit_error_total2_SSKM + length(bit_error_numb2_SSKM);
                
    end

    % BER for SC and MRC
    BER_SC(i) = bit_error_total_SC / (Data_Size * Iteration);
    BER_MRC(i) = bit_error_total_MRC / (Data_Size * Iteration);

    % BER for SSKM_SC and SSKM_MRC
    BER1_SSKM(i) = bit_error_total1_SSKM / (Data_Size * Iteration);
    BER2_SSKM(i) = bit_error_total2_SSKM / (Data_Size * Iteration);

    bit_error_total_SC = 0;
    bit_error_total_MRC = 0;
    bit_error_total1_SSKM = 0;
    bit_error_total2_SSKM = 0;
end


% Plot results
figure;
semilogy(SNR, BER_SC, '*-', SNR, BER_MRC, 'o-', SNR, BER1_SSKM, 's-', SNR, BER2_SSKM, 'diamond');
xlabel('SNR (dB)');
ylabel('BER');
legend('SC', 'MRC','SSKM_SC', 'SSKM_MRC');


function deSSK = deSSK(Data)
    deSSK = zeros(1, length(Data)/2);

    for i = 1 : length(Data)/2
        if Data(2*i-1) == 1
            deSSK(1, i) = 1;
        else
         deSSK(1, i) = 0;
        end
    end
end