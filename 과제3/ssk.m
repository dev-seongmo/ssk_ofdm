%% ssk funtcion

function modulated_signal = ssk_modulation(data)
    % SSK 변조
    modulated_signal = zeros(1, length(data)*2);
    for i = 1:length(data)
        if data(i) == 1
            modulated_signal(2*i-1) = 1; % 1인 경우 1로 설정
        else
            modulated_signal(2*i) = 1; % 0인 경우 1로 설정
        end
    end
end

function received_data = ssk_demodulation(received_signal)
    % SSK 복조
    received_data = zeros(1, length(received_signal)/2);
    for i = 1:length(received_data)
        if received_signal(2*i-1) == 1 && received_signal(2*i) == 0
            received_data(i) = 1;
        elseif received_signal(2*i-1) == 0 && received_signal(2*i) == 1
            received_data(i) = 0;
        elseif received_signal(2*i-1) == 0 && received_signal(2*i) == 0
            received_data(i) = 0;
        else
            received_data(i) = 1;
        end
    end
end