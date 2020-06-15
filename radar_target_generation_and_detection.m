clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_range = 200;
range_resolution = 1;

speed_of_light = 3e8;
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains constant
target_position = 50;
target_velocity = 0.00000;
 


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
% bandwidth = speed_of_light / (2 * range_resolution)
B = speed_of_light / (2 * range_resolution);
% Tchirp = 5.5 * 2 * range_max / speed_of_light
Tchirp = 5.5 * 2 * max_range / speed_of_light;
% slope = bandwidth / Tchirp
slope = B / Tchirp;


%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    target_position = target_position + target_velocity;
    r_t(i) = target_position;
    td(i) = r_t(i) * 2 / speed_of_light;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2 * pi * (fc * t(i) + slope * t(i)^2 / 2));
    receive_time = t(i) - td(i);
    
    Rx(i) = cos(2 * pi * (fc * receive_time + slope * receive_time^2 / 2));
    %doppler_ratio = 1 + target_velocity / speed_of_light
    %Rx(i) = cos(doppler_ratio * 2 * pi * (fc * receive_time + slope * receive_time^2 / 2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i) * Rx(i);
    %Mix(i) = cos(2 * pi * (2 * slope * r_t(i) * t(i) / speed_of_light + 2 * fc * target_velocity * t(i) / speed_of_light));
    
end

%figure ('Name','Range from First FFT')
%subplot(2,1,1);
%   plot(Tx);
%subplot(2,1,2);
%   plot(Rx);

%figure ('Name','Range from First FFT')
%subplot(2,1,1);
%   plot(Tx);
%subplot(2,1,2);
%   plot(Rx);
%   return
%   return
%% RANGE MEASUREMENT

 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix2d=reshape(Mix,[Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
range_fft = fft(Mix2d, Nr, 1);
normed_range_fft = range_fft / Nr;

 % *%TODO* :
% Take the absolute value of FFT output
abs_normed_range_fft = abs(normed_range_fft);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
pos_abs_normed_range_fft = abs_normed_range_fft(1:Nr/2 + 1, :);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1);

%for i=1:Nd
%    plot(pos_abs_normed_range_fft(:, i))
%    axis ([0 200 0 0.5]);
%    w = waitforbuttonpress;
%end

 % *%TODO* :
 % plot FFT output
first_fft = pos_abs_normed_range_fft(:, 1);
plot(first_fft)
axis ([0 200 0 0.5]);


%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix2d,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(Nr/2+1:Nr,1:Nd);
sig_fft2 = fftshift(sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM);

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis, range_axis, RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map
size_x = size(RDM, 1);
size_y = size(RDM, 2);

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
T = 6; % number_of_training_cells

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
G = 2; % number_of_guard_cells

% *%TODO* :
% offset the threshold by SNR value in dB
threshold = 3;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(size_x, size_y);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.
cfar = zeros(size_x, size_y);
number_of_training_cells = (1 + G + T)^2 - (1 + G)^2;

for x = (G + T + 1):(size_x - G - T)
    for y = (G + T + 1):(size_y - G - T)

        total_noise = 0;
        % 2. - 5. Determine the noise threshold by measuring it within the training cells
        for i = (-G - T):(G + T)
            for j = (-G - T):(G + T)
                if (i<-G || i>G) || (j<-G || j>G)
                    total_noise = total_noise + db2pow(RDM(x+i, y+j));
                end
            end
        end
        
        average_noise = total_noise / number_of_training_cells;
        % Convert noise level back to log scale
        noise_level(x, y) = pow2db(average_noise);

        % 6. Measuring the signal within the CUT
        cell_under_test_value = RDM(x, y);

        % 8. Filter the signal above the threshold
        if cell_under_test_value > (noise_level(x, y) + threshold)
            cfar(x, y) = 1;
        else
            cfar(x, y) = 0;
        end
    end
end





% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 








% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis, range_axis, cfar);
colorbar;


%figure,surf(doppler_axis, range_axis, noise_level);
%colorbar;



 
 