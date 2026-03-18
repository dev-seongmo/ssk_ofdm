clear clc;

% Main script
FFT_Size = 128;
GI_Size = FFT_Size / 4;
Modulation_Order = 4;
 
Iterations = 10000;
Multi_Path = 7;
Data_Size = FFT_Size*Modulation_Order;
BER = zeros(1,11);
Error_count = zeros(1,11);
for SNR  = 0:3:30

        for i = 1:Iterations
            
            data = randi([0, 1], 1, Data_Size);  
            
            
            modulated_data = base_mod(data, Modulation_Order);
            
            ofdm_data = ifft(modulated_data, FFT_Size)*sqrt(FFT_Size);
            cp_data = [ofdm_data(end-GI_Size+1:end),ofdm_data];
            h = Rayleigh_channel(Multi_Path);
            hx = conv(h,cp_data);
            Y = awgn_noise(hx,SNR);

            %수신
            cp_Y = Y(GI_Size+1:end-6);
            
            fft_data = fft(cp_Y,FFT_Size)/sqrt(FFT_Size);
            H = fft(h,FFT_Size);
            x = fft_data./H;
            
            data2 = base_demod(x,Modulation_Order);
            for k=1:Data_Size;
                if (data(k) ~= data2(k))
                Error_count((SNR/3)+1)=Error_count((SNR/3)+1)+1;
                end
            end
        end

        
        
end
BER = Error_count./(Data_Size*Iterations);


xSNR = 0:3:30;
semilogy(xSNR,BER,"-o");
title("BER performance");
legend("SISO-OFDM");
ylim([0.001 1]);
ylabel("BER");
xlabel("SNR(db)");
grid on;

