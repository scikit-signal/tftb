function [scalo,f,T,a,wt,wavescaled] = contwtgn(x,fmin,fmax,N,wave);

% CONTWTGN      [scalo,f,T,a,wt,wavescaled] = contwtgn(x,fmin,fmax,N,wave) 
%				computes a continuous wavelet transform.
%
% Input:
%	x		signal (in time) to be analyzed          
%	[fmin,fmax]  	respectively lower and upper frequency bounds of
%           		the analysis (in cycles/sec).
%	[N] 		number of analyzed voices
%	[wave] 		specifies the analyzing wavelet
%                       An order "wave" derivative of the Gaussian is chosen
%
% Output:    
%    scalo 		scalogram (squared magnitude of WT)
%	 f 		frequency samples (geometrically sampled between FMAX 
%			and FMIN).
%        T              time samples
%	 a 		scale vector (geometrically sampled between 
%			1 and FMAX/FMIN)
%	 wt 		coefficient of the wavelet transform. 
%	        	X-axis corresponds to time (uniformly sampled), Y-axis
%          	 	corresponds to frequency (or scale) samples  
%	       		(geometrically sampled between Fmin (resp. Fmax/Fmin) 
%			and Fmax (resp. 1)
%           		First row of WT corresponds to the lowest analyzed 
%			frequency.
% wavescaled 		when the analyzing wavelet is Morlet or Mexican hat, 
%			wavescaled = wave. For an aritrary band-pass analyzing 
%			function, wavescaled contains columnwise the (N) 
%			scaled version of it
%
% Example:    	S = altes(256,0.1,0.45,10000) ;
%             	[scalo,f,T,a,wt,wavescaled] = contwtgn(S,0.01,0.5,128,8) ;
%
% See also:
%

% Permission is granted for use and non-profit distribution providing that this
% notice be clearly maintained. The right to distribute any portion for profit
% or as part of any commercial product is specifically reserved for the author.


% CHECK INPUT FORMATS

[xr,xc] = size(x) ;
if xr ~= 1 & xc ~= 1
  error('1-D signals only')
elseif xc == 1
  x = x.' ;
end

T = 1 : length(x) ;

% DEFAULT VALUES

nt = size(x,2) ;
if exist('wave') == 0 
  wave = 2 ;
end

if nargin == 1
  XTF = fft(fftshift(x)) ;
  sp = (abs(XTF(1:nt/2))).^2 ;
  f = linspace(0,0.5,nt/2+1) ; f = f(1:nt/2) ;
  plot(f,sp) ; grid ;
  xlabel('frequency');
  title('Analyzed Signal Spectrum') ;
  fmin = input('lower frequency bound = ') ;
  fmax = input('upper frequency bound = ') ;
  N = input('Frequency samples = ') ;
  fmin_s = num2str(fmin) ; fmax_s = num2str(fmax) ; 
  N_s = num2str(N) ;
  disp(['frequency runs from ',fmin_s,' to ',fmax_s,' over ',N_s,' points']) ;
end
if nargin == 4 | nargin == 5
  if fmin >= fmax
    error('fmax must be greater than fmin') ;
  end
end

f = logspace(log10(fmax),log10(fmin),N) ;
a = logspace(log10(1),log10(fmax/fmin),N) ; amax = max(a) ;

for ptr = 1:N
  ha = gaussn(f(ptr),wave) ; nha = (length(ha)-1)/2 ;
  detail = conv(x,ha) ;
  wt(ptr,:) = detail(nha+1:nha+nt) ;
end  

wavescaled = wave ;


%%%% pour etre compatible avec le format de donnees de TFTB %%%%%

wt = flipud(wt) ;
a = flipud(a(:)) ;
f = flipud(f(:)) ;
scalo = (wt.*conj(wt)) ;


%%%% pour etre compatible avec le format de donnees de TFTB %%%%%


%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
