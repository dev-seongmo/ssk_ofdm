% Main script
FFT_Size = 128;
GI_Size = FFT_Size / 4;
Modulation_Orders = 4;
SNRs = 0:3:30;  
Iterations = 10000;
Multi_Path = 7;

BER = zeros(length(SNRs), length(Modulation_Orders));

for m = 1:length(Modulation_Orders)
    Modulation_Order = Modulation_Orders(m);
    Data_Size = FFT_Size * Modulation_Order;

    for s = 1:length(SNRs)
        SNR = SNRs(s);
        error_count = 0;
        total_bits = 0;

        for iter = 1:Iterations
            % 데이터 생성
            data = randi([0 Modulation_Order^2-1], 1, Data_Size);

            % 변조
            modulated_data = base_mod(data, Modulation_Order,Data_Size);

            % IFFT 수행 (OFDM 변환)
            ofdm_data = ifft(modulated_data, FFT_Size);

            % 가드 인터벌 추가
            tx_signal = [ofdm_data(end-GI_Size+1:end) ofdm_data];

            % 채널 통과
            h = Rayleigh_channel(Multi_Path);
            rx_signal = conv(h, tx_signal);
            rx_signal = rx_signal(1:length(tx_signal)); % 트림

            % 잡음 추가
            rx_signal = awgn_noise(rx_signal, SNR);

            % 가드 인터벌 제거
            rx_signal = rx_signal(GI_Size+1:end);

            % FFT 수행 (OFDM 변환)snr
            received_data = fft(rx_signal, FFT_Size);

            % 복조
            demodulated_data = base_demod(received_data, Modulation_Order,Data_Size);

            % 오류 계산
            error_count = error_count + sum(data ~= demodulated_data);
            total_bits = total_bits + Data_Size;
        end

        BER(s, m) = error_count / total_bits;
    end
end

% 결과 시각화
figure;
for m = 1:length(Modulation_Orders)
    semilogy(SNRs, BER(:, m), 'DisplayName', ['Modulation Order: ' num2str(Modulation_Orders(m))]);
    hold on;
end
xlabel('SNR (dB)');
ylabel('BER');
legend show;
title('BER vs SNR for different Modulation Orders');
grid on;
