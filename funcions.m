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
	fmax = 1024;

	if N > length(x) - k;
		N = length(x) - k;
		fprintf('Time is soo huge!. \n')
	end

	%aby se dodrzel vzorkovaci teorem musi byt N > 2*max(freq)
	%fmax -> 524 Mhz
	if N < fmax %add freq max
		n = fmax - N;
	end

	ff=0:fs/N:fs-fs/N;


	figure(1);
	subplot(511);
	plot(x);

	X = fft(x(1+k:N+k), N + n);
	Xabs = abs(X);
	subplot(512);
	plot(x(1+k:N+k));
	subplot(513);
	plot(ff, Xabs);

	d = func.maxval(Xabs, 0.8);
	d = ff(d)

end

function [shift] = find_start(x, fs, t)
	func =  funcions;
	shift = func.maxval(x, 0.9);
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
	a = 2^(1/12);
	freq = [440/a^9, 440/a^8, 440/a^7, 440/a^6, 440/a^5, 440/a^4, 440/a^3, 440/a^2, 440/a, 440, 440*a, 440*a^2, 440*a^3];
	[c i] = min(abs(freq-x));
	d = freq(i);
	tone = tones(i,:);
end
