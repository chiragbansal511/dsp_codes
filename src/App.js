import React from 'react';
import './App.css';

const codes = {
  '1bnd': `% Prompt the user for frequency and duty cycle
freq = input('Enter the frequency of the square wave (Hz): ');
duty_cycle = input('Enter the duty cycle of the square wave (0-100%): ');

fs = 1000;
t_end = 1;
t = 0:1/fs:t_end-1/fs;

x_cont = square(2 * pi * freq * t, duty_cycle);
x_discrete = square(2 * pi * freq * t, duty_cycle);

figure;
subplot(2, 1, 1);
plot(t, x_cont);
title('Continuous-Time Square Wave');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(2, 1, 2);
stem(t, x_discrete, 'filled');
title('Discrete-Time Square Wave');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;`,

  '1 c': `% Prompt user for the number of sinusoids
num_sinusoids = input('Enter the number of sinusoids: ');

frequencies = zeros(1, num_sinusoids);
amplitudes = zeros(1, num_sinusoids);
phases = zeros(1, num_sinusoids);

for i = 1:num_sinusoids
    fprintf('Sinusoid %d:\n', i);
    frequencies(i) = input('  Frequency (Hz): ');
    amplitudes(i) = input('  Amplitude: ');
    phases(i) = input('  Phase (radians): ');
end

fs = 1000;
t_end = 1;
t = 0:1/fs:t_end-1/fs;

x = zeros(size(t));

for i = 1:num_sinusoids
    x = x + amplitudes(i) * sin(2 * pi * frequencies(i) * t + phases(i));
end

figure;
plot(t, x);
title('Composite Sinusoidal Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;`,

  '2': `% Define two discrete sequences
x = [1, 2, 3, 4];
h = [1, 0, -1];

M = length(x);
N = length(h);

L = M + N - 1;

y = zeros(1, L);

for n = 1:L
    convolution_sum = 0;
    for k = 1:M
        if (n - k + 1) > 0 && (n - k + 1) <= N
            convolution_sum = convolution_sum + x(k) * h(n - k + 1);
        end
    end
    y(n) = convolution_sum;
end

disp('Convolution result:');
disp(y);

figure;
subplot(3, 1, 1);
stem(x);
title('Sequence x[n]');
xlabel('n');
ylabel('x[n]');
grid on;

subplot(3, 1, 2);
stem(h);
title('Sequence h[n]');
xlabel('n');
ylabel('h[n]');
grid on;

subplot(3, 1, 3);
stem(y);
title('Convolution y[n]');
xlabel('n');
ylabel('y[n]');
grid on;`,

  '3 a': `% Input sequence from the user
x = input('Enter the input sequence: ');

N = length(x);

X = zeros(1, N);

for k = 0:N-1
    for n = 0:N-1
        X(k+1) = X(k+1) + x(n+1) * exp(-1j * 2 * pi * k * n / N);
    end
end

magnitude_spectrum = abs(X);
phase_spectrum = angle(X);

figure;
subplot(2, 1, 1);
stem(0:N-1, magnitude_spectrum, 'filled');
title('Magnitude Spectrum');
xlabel('Frequency Bin');
ylabel('|X[k]|');
grid on;

subplot(2, 1, 2);
stem(0:N-1, phase_spectrum, 'filled');
title('Phase Spectrum');
xlabel('Frequency Bin');
ylabel('Phase of X[k] (radians)');
grid on;`,

  '3 b': `% Input sequence from the user
x = input('Enter the input sequence: ');

N = length(x);

X = zeros(1, N);

for k = 0:N-1
    for n = 0:N-1
        X(k+1) = X(k+1) + x(n+1) * exp(-1j * 2 * pi * k * n / N);
    end
end

x_idft = zeros(1, N);
for n = 0:N-1
    for k = 0:N-1
        x_idft(n+1) = x_idft(n+1) + X(k+1) * exp(1j * 2 * pi * k * n / N);
    end
    x_idft(n+1) = x_idft(n+1) / N;
end

disp('Original sequence:');
disp(x);

disp('Reconstructed sequence from IDFT:');
disp(real(x_idft));

figure;
subplot(2, 1, 1);
stem(0:N-1, x, 'filled');
title('Original Sequence');
xlabel('n');
ylabel('x[n]');
grid on;

subplot(2, 1, 2);
stem(0:N-1, real(x_idft), 'filled');
title('Reconstructed Sequence from IDFT');
xlabel('n');
ylabel('x_{IDFT}[n]');
grid on;`,

  '4 a': `% Input sequences from the user
x = input('Enter the first input sequence x[n]: ');
h = input('Enter the second input sequence h[n]: ');

N = max(length(x), length(h));

x = [x, zeros(1, N - length(x))];
h = [h, zeros(1, N - length(h))];

y_circular = zeros(1, N);

for n = 0:N-1
    for k = 0:N-1
        y_circular(n+1) = y_circular(n+1) + x(k+1) * h(mod(n - k, N) + 1);
    end
end

disp('Circular convolution result (manual method):');
disp(y_circular);`,

  '4 b': `% Define your circular convolution result and the length of the original sequences
circular_conv = [1 2 3 4 5 6 7 8];
N = 4;

L = length(circular_conv);

padding_length = 2 * N;

circular_conv_padded = [circular_conv, zeros(1, padding_length - L)];

linear_conv = zeros(1, padding_length);

for n = 1:padding_length
    for k = 1:padding_length
        linear_conv(n) = linear_conv(n) + circular_conv_padded(mod(n - k, L) + 1) * circular_conv_padded(mod(k - n, L) + 1);
    end
end

disp('Linear Convolution Result:');
disp(linear_conv);`,

  '4 c': `% Define your sequences
x = [1 2 3 4];
h = [1 1 1 1];

Nx = length(x);
Nh = length(h);

N = Nx + Nh - 1;

x_padded = [x, zeros(1, N - Nx)];
h_padded = [h, zeros(1, N - Nh)];

linear_conv = zeros(1, N);

for n = 1:N
    for k = 1:N
        if (n - k + 1) > 0
            linear_conv(n) = linear_conv(n) + x_padded(k) * h_padded(n - k + 1);
        end
    end
end

circular_conv = linear_conv(1:Nx);

disp('Circular Convolution Result:');
disp(circular_conv);`,

  '4 d': `% Define your sequences
x = [1 2 3 4]; % Example sequence x
h = [1 1 1 1]; % Example sequence h

% Length of sequences
Nx = length(x);
Nh = length(h);

% Determine the length for zero-padding
N = max(Nx, Nh);

% Zero-pad the sequences to the length of the maximum of Nx and Nh
x_padded = [x, zeros(1, N - Nx)];
h_padded = [h, zeros(1, N - Nh)];

% Define the DFT function
function X = myDFT(x)
    N = length(x);
    X = zeros(1, N);
    for k = 0:N-1
        sum = 0;
        for n = 0:N-1
            angle = 2 * pi * k * n / N;
            sum = sum + x(n+1) * (cos(angle) - 1i * sin(angle));
        end
        X(k+1) = sum;
    end
end

% Define the IDFT function
function x = myIDFT(X)
    N = length(X);
    x = zeros(1, N);
    for n = 0:N-1
        sum = 0;
        for k = 0:N-1
            angle = 2 * pi * k * n / N;
            sum = sum + X(k+1) * (cos(angle) + 1i * sin(angle));
        end
        x(n+1) = sum / N;
    end
end

% Compute the DFT of each sequence
X = myDFT(x_padded);
H = myDFT(h_padded);

% Multiply the DFTs element-wise
Y = X .* H;

% Compute the IDFT of the product to get the circular convolution result
y = myIDFT(Y);

% Display the result
disp('Circular Convolution Result:');
disp(real(y)); % Use real() to discard small imaginary parts due to numerical precision`,

  '5': `% Define the window length
N = 64; % Length of the window

% Define the fixed windows
rectangular_window = ones(1, N); % Rectangular window
hamming_window = hamming(N)'; % Hamming window
hanning_window = hanning(N)'; % Hanning window

% Define the DFT function
function X = myDFT(x)
    N = length(x);
    X = zeros(1, N);
    for k = 0:N-1
        sum = 0;
        for n = 0:N-1
            angle = 2 * pi * k * n / N;
            sum = sum + x(n+1) * (cos(angle) - 1i * sin(angle));
        end
        X(k+1) = sum;
    end
end

% Define the IDFT function
function x = myIDFT(X)
    N = length(X);
    x = zeros(1, N);
    for n = 0:N-1
        sum = 0;
        for k = 0:N-1
            angle = 2 * pi * k * n / N;
            sum = sum + X(k+1) * (cos(angle) + 1i * sin(angle));
        end
        x(n+1) = sum / N;
    end
end

% Compute the frequency domain responses
rectangular_freq = abs(myDFT(rectangular_window));
hamming_freq = abs(myDFT(hamming_window));
hanning_freq = abs(myDFT(hanning_window));

% Plot the time domain response of the windows
figure;
subplot(3, 1, 1);
plot(rectangular_window);
title('Rectangular Window (Time Domain)');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(hamming_window);
title('Hamming Window (Time Domain)');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(hanning_window);
title('Hanning Window (Time Domain)');
xlabel('Sample Index');
ylabel('Amplitude');

% Plot the frequency domain response of the windows
figure;
subplot(3, 1, 1);
plot(rectangular_freq);
title('Rectangular Window (Frequency Domain)');
xlabel('Frequency Bin');
ylabel('Magnitude');

subplot(3, 1, 2);
plot(hamming_freq);
title('Hamming Window (Frequency Domain)');
xlabel('Frequency Bin');
ylabel('Magnitude');

subplot(3, 1, 3);
plot(hanning_freq);
title('Hanning Window (Frequency Domain)');
xlabel('Frequency Bin');
ylabel('Magnitude');`
};

function App() {
  const copyToClipboard = (text) => {
    navigator.clipboard.writeText(text).then(() => {
      // alert('Code copied to clipboard!');
    }).catch(err => {
      console.error('Failed to copy: ', err);
    });
  };

  return (
    <div className="App">
      <h1>MATLAB Code Snippets</h1>
      <div className="button-container">
        {Object.keys(codes).map((key) => (
          <button className="code-button" key={key} onClick={() => copyToClipboard(codes[key])}>
            {key}
          </button>
        ))}
      </div>
    </div>
  );
}

export default App;

