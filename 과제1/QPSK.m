% 64개의 랜덤한 0과 1로 구성된 배열 생성
datasize = 5096 ;  
order =4;
data = randi([0, 1], 1, datasize);  
streamSym = zeros(1, datasize / order);
H = (randn(1) + 1j * randn(1)) / sqrt(2);
N = (randn(1,datasize/order) + 1j * randn(1,datasize/order)) / sqrt(2);
for idx = 1:datasize/order
    if (data(4*idx-3))
        a=-1;
    else
        a=+1;
    end
    if (data(4*idx-2))
        b=1;
    else
        b=-1;
    end
    if (data(4*idx-1))
        a2=3;
    else
        a2=1;
    end
    if (data(4*idx))
        b2=3;
    else
        b2=1;
    end
    streamSym(idx)=a*a2+b*b2*1j;

end


 Y= H*streamSym + N;
  scatter(real(Y), imag(Y))
   hold on

    scatter(real(streamSym), imag(streamSym),"blue","filled");
 hold on

    
