%--------------------------------------------------------------------------
% General parameters for the simulation
%--------------------------------------------------------------------------
T = 1;                 % Symbol duration
N = 16;                % Number of carriers
M = 16;                % Constellation order (16-QAM)
m = log2(M);           % Bits per symbol
nBlocks = 1e5;         % Number of OFDM blocks
nSymb = nBlocks * N;   % Total number of symbols
nBits = nSymb * m;     % Total number of bits
tAssign = 'gray';      % Type of binary assignment ('gray', 'bin')

a = 9/10;
d = [a^0 a^1 a^2 a^3 a^4 a^5];  % Channel impulse response



% Noise variance
varNoise = 4; 

%cyclic prefix lengths to test
C_values = 1:10;  

%--------------------------------------------------------------------------
% Transmission and Error Probability Calculation
%--------------------------------------------------------------------------
Pe_mean = zeros(1, length(C_values)); 

for idx = 1:length(C_values)
    C = C_values(idx); % Current cyclic prefix length
    
    % Generación de bits (Generation of Bits) 
    B = randi([0 1], nBits, 1);          
    % Codificación de bits en símbolos (Symbols encoded from bits)
    A = qammod(B, M, tAssign, 'InputType', 'bit'); 
    
    % Conversión serie a paralelo / Serial to parallel conversion
    Ak = reshape(A, N, nBlocks);
    sb = N / sqrt(T) * ifft(Ak, N, 1);

    %add cyclic prefix (C)
    sb_cp = [sb(end-C+1:end, :); sb];

    % Conversión paralelo a serie / Parallel to serial conversion
    s = reshape(sb_cp, 1, []);

    %--------------------------------------------------------------------------
    % Transmission
    %--------------------------------------------------------------------------
    Kd = length(d) - 1; 
    v = conv([zeros(1, Kd), s], d);
    v = v(Kd+1:Kd+length(s));
    v = v + sqrt(varNoise / 2) * (randn(size(v)) + 1j * randn(size(v)));

    %-------------------------------------------------
    % Demodulation
    %--------------------------------------------------
    v_blocks = reshape(v, N + C, nBlocks);
    v_no_cp = v_blocks(C+1:end, :);
    qk = 1 / sqrt(T) * fft(v_no_cp, N, 1);

    %--------------------------------------------------------------------------
    % Visualization
    %--------------------------------------------------------------------------
    p = OFDM(d, N, C); 
    Pe = zeros(1, N);

    for k = 1:N
        pk = p{k, k};
        qn = qk(k, :) / pk(1); 
        Bke = qamdemod(qn, M, tAssign, 'OutputType', 'bit'); 
        Ake = qammod(Bke, M, tAssign, 'InputType', 'bit'); 
        if C == 5
            Pe(k)=length(find(Ake~=Ak(k,:)))/length(Ake);
            disp(['  Pe (k=',num2str(k-1),') = ',num2str(Pe(k))])
        end
        Pe(k) = sum(Ake ~= Ak(k, :)) / length(Ake);  
    end

    Pe_mean(idx) = mean(Pe);
    
    disp(['Cyclic Prefix Length C = ', num2str(C), ', Mean Pe = ', num2str(Pe_mean(idx))]);
end

%--------------------------------------------------------------------------
% Plotting results
%--------------------------------------------------------------------------
figure;
plot(C_values, Pe_mean, '-o');
grid on;
xlabel('Cyclic Prefix Length (C)');
ylabel('Mean Symbol Error Probability (Pe)');
title('Mean Pe vs Cyclic Prefix Length');