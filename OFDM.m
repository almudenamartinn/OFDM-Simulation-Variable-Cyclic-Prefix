function p = OFDM(d, N, C)
    % OFDM_p generates a processing matrix based on the channel impulse response
    % d: Channel impulse response (length Kd+1)
    % N: Number of subcarriers
    % C: Length of cyclic prefix

    % Length of the channel impulse response
    Kd = length(d) - 1; 

    % Frequency response of the channel
    H = fft(d, N);

    % Create a cell array to store the channel frequency responses
    p = cell(N, N);

    % Loop to populate the processing matrix p
    for i = 1:N
        for j = 1:N
            % Here we use the channel frequency response H and apply it
            % to each element of the matrix p based on subcarrier indices
            p{i, j} = H(i) * H(j);
        end
    end
end
