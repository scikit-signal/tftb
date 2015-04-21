function [y,iflaw]=fmsin(N,fnormin,fnormax,period,t0,fnorm0,pm1);
%FMSIN	Signal with sinusoidal frequency modulation.
%	[Y,IFLAW]=FMSIN(N,FNORMIN,FNORMAX,PERIOD,T0,FNORM0,PM1)
%	generates a frequency modulation with a sinusoidal frequency.
%	This sinusoidal modulation is designed such that the instantaneous
%	frequency at time T0 is equal to FNORM0, and the ambiguity 
%	between increasing or decreasing frequency is solved by PM1.
%
%	N       : number of points.
%	FNORMIN : smallest normalized frequency          (default: 0.05) 
%	FNORMAX : highest normalized frequency           (default: 0.45)
%	PERIOD  : period of the sinusoidal fm            (default: N   )
%	T0      : time reference for the phase           (default: N/2 )
%	FNORM0  : normalized frequency at time T0        (default: 0.25)
%	PM1     : frequency direction at T0 (-1 or +1)	 (default: +1  )
%	Y       : signal
%	IFLAW   : its instantaneous frequency law (optional).
%
%	Example: 
%	 z=fmsin(140,0.05,0.45,100,20,0.3,-1.0);plot(real(z));
%
%	See also FMCONST, FMLIN, FMODANY, FMHYP, FMPAR, FMPOWER.

%	F. Auger, July 1995.
%	Copyright (c) 1996 by CNRS (France).
%
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

if (nargin == 0),
 error ( 'At least 1 parameter required' ) ;
elseif (nargin == 1),
 fnormin = 0.05; fnormax = 0.45; period = N; 
 t0= round(N/2); fnorm0 = 0.25; pm1=+1;
elseif (nargin == 2),
 fnormax = 0.45; period = N; 
 t0 = round(N/2); fnorm0 = 0.5*(fnormin+fnormax); pm1=+1;
elseif (nargin == 3),
 period = N; t0 = round(N/2); 
 fnorm0 = 0.5*(fnormin+fnormax); pm1=+1;
elseif (nargin == 4),
 t0 = round(N/2); fnorm0 = 0.5*(fnormin+fnormax); pm1=+1;
elseif (nargin == 5),
 fnorm0 = 0.5*(fnormin+fnormax); pm1=+1;
elseif (nargin == 6),
 pm1= +1;
elseif (nargin==7),
 if (abs(pm1)~=1),
  error('pm1 must be equal to -1 or +1');
 end;
end;

if (N <= 0),
 error ('The signal length N must be strictly positive' );
elseif (abs(fnormin) > 0.5)|(abs(fnormax) > 0.5)|(abs(fnorm0) > 0.5),
 error ('fnormin, fnormax and fnorm0 must be between -0.5 and 0.5') ;
elseif (fnormin > fnormax)
 error ('fnormin must be lower than fnormax');
elseif (fnormin > fnorm0)|(fnorm0 > fnormax),
 error ('fnorm0 must be between fnormin and fnormax') ;
else
 fnormid=0.5*(fnormax+fnormin);
 delta  =0.5*(fnormax-fnormin);
 phi=-pm1*acos((fnorm0-fnormid)/delta);
 time=(1:N)-t0;
 phase=2*pi*fnormid*time+delta*period*(sin(2*pi*time/period+phi)-sin(phi));
 y = exp(j*phase).';
 if (nargout==2)
  iflaw=fnormid+delta*cos(2*pi*time'/period+phi);
 end
end
