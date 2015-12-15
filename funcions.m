%%
function func=funcions
  func.dist = @dist;
  func.fft = @fft_w;
  func.rec = @recognize;
  func.find = @find_start;
  func.maxval = @max_val;
  func.map = @map;
end

function [d] = dist(a,x)
    d = a-x;
end

function [d] = fft_w(x, fs, t, k)
	func =  funcions;
	n = 0;
	N = fs*t; %delka segmenty pro delky ms
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

	x = filter(b, 1, x);

	ff=0:fs/(N+n):fs-fs/(N+n);

	figure(1);
	subplot(311);
	plot(x);

	wN=hamming(N);
	swN=x(1+k:N+k).*wN(1:N);
	X = fft(swN, N);
	Xabs = abs(X);
	subplot(312);
	plot(x(1+k:N+k));
	subplot(313);
	limit = floor(length(ff)/50);
	ff = ff(1:limit);
	Xabs = Xabs(1:limit);
	plot(ff, Xabs);

	d = func.maxval(Xabs, 0.7);
	d = ff(d);
end

function [shift] = find_start(x, fs, t)
	func =  funcions;
	shift = func.maxval(x);
end

function [shift] = max_val(x, k)
	if nargin == 1
		k = 1;
	end
	shift = find(x >= k*max(x));
	shift = shift(1);
end

function [d] = recognize(x, fs, t, k)
	func =  funcions;
	k = func.find(x, fs, t);
	d = func.fft(x, fs, t, k);
	d = func.map(d);
end


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
