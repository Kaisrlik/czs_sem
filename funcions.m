%% functions
function func=funcions
  func.fft = @fft_w;
  func.rec = @recognize;
  func.find = @find_start;
  func.maxval = @max_val;
  func.map = @map;
  func.mach = @mach;
end

function [freq] = fft_w(x, fs, t, k)
	func =  funcions;
	n = 0;
	N = fs*t;
	fmax = 2048;

	if N > length(x) - k;
		N = length(x) - k;
		fprintf('Time is soo huge!. \n')
	end

	%aby se dodrzel vzorkovaci teorem musi byt N > 2*max(freq)
	%fmax -> 524 Mhz
	if N < fmax %add freq max
		n = fmax - N;
	end

	wn = [200/(fs/2), 600/(fs/2)];
	b = fir1(64, wn, 'bandpass');
%	figure(5);
%	zplane(b);
%	figure(3);
%	freqz(b, 1);
%	figure(4);
%	impz(b, 1);
	x = filter(b, 1, x); % XXX useless

	ff=0:fs/(N+n):fs-fs/(N+n);

	figure(1);
	subplot(411);
	plot(x);
	title('Whole signal');
	ylabel('Ampl');
	xlabel('Time');

	wN=hamming(N);
	swN=x(1+k:N+k).*wN(1:N);
	X = fft(swN, N+n);
	Xabs = abs(X);
%	figure(2);
	subplot(412);
	plot(x(1+k:N+k));
	title('The interesting part of signal');
	ylabel('Ampl');
	xlabel('Time');
%	figure(3);
	subplot(413);
	plot(ff, Xabs);
	title('FFT of signal above');
	ylabel('Magnitude');
	xlabel('Frequency');
%	figure(4);
	subplot(414);
	limit = floor(length(ff)/50);
	ff = ff(1:limit);
	Xabs = Xabs(1:limit);
	plot(ff, Xabs, '-o');
	title('Zoom in FFT');
	ylabel('Magnitude');
	xlabel('Frequency');

	freq = func.mach(Xabs, ff);
end

% finding of main 
function [freq] = mach(X, ff)
	func =  funcions;
	freq = -1;
	[pks, locs] = findpeaks(X);
	%locs = locs(find(pks > 0.1*max(pks))) % add some constrains
	pksf = ff(locs);
	A = ff(func.maxval(X, 0.55));

	cont = ismember(A, pksf);

	if length(A) > 1
		for i = 1 : length(A);
			if cont(i) == 0
				continue;
			end
			[c] = min(abs(pksf-A(i)));
			if c > 5
				continue;
			end
			[c] = min(abs(pksf-A(i)*2));
			if c < 20 % XXX this constrain can be smaller?
				freq = A(i);
				return;
			end
		end
	end

	fprintf('No mach!. \n');
	freq = A(find(cont == 1));
	freq = freq(1);
end

%
function [shift] = find_start(x, fs, t)
	func =  funcions;
	shift = func.maxval(x);
	shift = shift(1);
end

% finding maximum passing limitation of k
function [shift] = max_val(x, k)
	if nargin == 1
		k = 1;
	end
	shift = find(x >= k*max(x));
end

function [tone, freq] = recognize(x, fs, t, k)
	func =  funcions;
	k = func.find(x, fs, t);
	freq = func.fft(x, fs, t, k);
	tone = func.map(freq);
end

% mapping freq to nearest one
function [tone] = map(x)
	tones = char('c', 'cis', 'd', 'dis', 'e', 'f', 'fis', 'g', 'gis', 'a', 'ais', 'h', 'c*');
	s = 2^(1/12);
	a = 440;
	freq = [a/s^9, a/s^8, a/s^7, a/s^6, a/s^5, a/s^4, a/s^3, a/s^2, a/s, a, a*s,a*s^2, a*s^3];
	[c i] = min(abs(freq-x));
	d = freq(i);
	tone = tones(i,:);
	if c > 100
		tone = 'none';
	end
end
