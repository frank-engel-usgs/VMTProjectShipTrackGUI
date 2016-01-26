function ddvec = dms2dd(d,m,s,n)

%DMS2DD Converts a [deg min sec] matrix to a decimal degrees vector format
%
%  dms = DMS2DD(d,m,s) converts a deg:min:sec matrix into a vector
%  format.  The vector format is dd = deg + min/60 + sec/3600.
%  This allows d,m,s triple to be compressed into a single value,
%  which can then be employed similar to a degree or radian vector.
%  The inputs d, m and s must be of equal size.  Minutes and
%  second must be between 0 and 60.
%
%  dms = DMS2DD(mat) assumes and input matrix of [d m s].  This is
%  useful only for single column vectors for d, m and s.
%
%  dms = DMS2DD(d,m) and dms = MAT2DMS([d m]) assume that seconds
%  are zero, s = 0.
%
%  dms = DMS2DD(d,m,s,n) uses n as the accuracy of the seconds
%  calculation.  n = -2 uses accuracy in the hundredths position,
%  n = 0 uses accuracy in the units position.  Default is n = -5.
%  For further discussion of the input n, see ROUNDN.
%
%  See also DMS2MAT

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.8 $    $Date: 1998/08/10 17:47:54 $
%  Revised by: Frank L. Engel, USGS 2014/05/27


if nargin == 0
   error('Incorrect number of arguments')

elseif nargin==1
   if size(d,2)== 3
       s = d(:,3);   m = d(:,2);   d = d(:,1);
   elseif size(d,2)== 2
       m = d(:,2);   d = d(:,1);   s = zeros(size(d));
   elseif size(d,2) == 0
       d = [];   m = [];   s = [];
   else
       error('Single input matrices must be n-by-2 or n-by-3.');
   end
   n = -5;

elseif nargin == 2
   s = zeros(size(d));
   n = -5;

elseif nargin == 3
   n = -5;
end

%  Test for empty arguments

if isempty(d) & isempty(m) & isempty(s);  ddvec = [];  return;  end

%  Don't let seconds be rounded beyond the tens place.
%  If you did, then 55 seconds rounds to 100, which is not good.

if n == 2;  n = 1;   end

%  Complex argument tests

if any([~isreal(d) ~isreal(m) ~isreal(s)])
    warning('Imaginary parts of complex ANGLE argument ignored')
	d = real(d);   m = real(m);   s = real(s);
end

%  Dimension and value tests

if ~isequal(size(d),size(m),size(s))
    error('Inconsistent dimensions for input arguments')
elseif any(rem(d(~isnan(d)),1) ~= 0 | rem(m(~isnan(m)),1) ~= 0)
    error('Degrees and minutes must be integers')
end

if any(abs(m) > 60) | any (abs(m) < 0)       %  Actually algorithm allows for
    error('Minutes must be >= 0 and < 60')   %  up to exactly 60 seconds or
                                             %  60 minutes, but the error message
elseif any(abs(s) > 60) | any(abs(s) < 0)    %  doesn't reflect this so that angst
    error('Seconds must be >= 0 and < 60')   %  is minimized in the user docs
end

%  Ensure that only one negative sign is present and at the correct location

if any((s<0 & m<0) | (s<0 & d<0) | (m<0 & d<0) )
    error('Multiple negative entries in a DMS specification')
elseif any((s<0 & (m~=0 | d~= 0)) | (m<0 & d~=0))
    error('Incorrect negative DMS specification')
end

%  Construct a sign vector which has +1 when
%  angle >= 0 and -1 when angle < 0.  Note that the sign of the
%  angle is associated with the largest nonzero component of d:m:s

negvec = (d<0) | (m<0) | (s<0);
signvec = ~negvec - negvec;

%  Convert to all positive numbers.  Allows for easier
%  adjusting at 60 seconds and 60 minutes

d = abs(d);     m = abs(m);    s = abs(s);

%  Truncate seconds to a specified accuracy to eliminate round-off errors

[s,msg] = roundn(s,n);
if ~isempty(msg);   error(msg);   end

%  Adjust for 60 seconds or 60 minutes. If s > 60, this can only be
%  from round-off during roundn since s > 60 is already tested above.
%  This round-off effect has happened though.

indx = find(s >= 60);
if ~isempty(indx);   m(indx) = m(indx) + 1;   s(indx) = 0;   end

%  The user can not put minutes > 60 as input.  However, the line
%  above may create minutes > 60 (since the user can put in m == 60),
%  thus, the test below includes the greater than condition.

indx = find(m >= 60);
if ~isempty(indx);   d(indx) = d(indx) + 1;   m(indx) = m(indx)-60;   end

%  Construct the dms vector format

ddvec = signvec .* (d + m./60 + s./3600);

%%%%%%%%%%%%%%%%
% Subfunctions %
%%%%%%%%%%%%%%%%
function [x,msg] = roundn(x,n)

%ROUNDN  Rounds input data at specified power of 10
%
%  y = ROUNDN(x) rounds the input data x to the nearest hundredth.
%
%  y = ROUNDN(x,n) rounds the input data x at the specified power
%  of tens position.  For example, n = -2 rounds the input data to
%  the 10E-2 (hundredths) position.
%
%  [y,msg] = ROUNDN(...) returns the text of any error condition
%  encountered in the output variable msg.
%
%  See also ROUND

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:48:04 $

msg = [];   %  Initialize output

if nargin == 0
    error('Incorrect number of arguments')
elseif nargin == 1
    n = -2;
end

%  Test for scalar n

if max(size(n)) ~= 1
   msg = 'Scalar accuracy required';
   if nargout < 2;  error(msg);  end
   return
elseif ~isreal(n)
   warning('Imaginary part of complex N argument ignored')
   n = real(n);
end

%  Compute the exponential factors for rounding at specified
%  power of 10.  Ensure that n is an integer.

factors  = 10 ^ (fix(-n));

%  Set the significant digits for the input data

x = round(x * factors) / factors;