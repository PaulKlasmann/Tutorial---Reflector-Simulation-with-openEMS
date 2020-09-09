function h = plotFFcocx(nf2ff,varargin)
%  h = plotFFcocx(nf2ff,varargin)
%
%  plot co- and cross-pol far field patterns in dBi
%
% input:
%   nf2ff:      output of CalcNF2FF
%
% variable input:
%   'freq_index':  - use the given frequency index, see nf2ff.freq
%                  - default is 1
%   'xaxis':       - 'phi' (default) or 'theta'
%   'param':       - array positions of parametric plot
%                  - if xaxis='phi', theta is parameter, and vice versa
%                  - default is 1
%
%   example:
%       plotFFdB(nf2cocx, 'freq_index', 2, ...
%                       'xaxis', 'phi', 'param', [1 46 91])
%
%       see examples/NF2FF/infDipol.m
%
% See also CalcNF2FF, plotFF3D, polarFF, plotFFdB
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig, Stefan Mahr
% Based up plotFFdB and modified for co- and cross-pol by Bruce Veidt

% defaults
freq_index = 1;
xaxis = 'phi';
param = 1;

for n=1:2:numel(varargin)
    if (strcmp(varargin{n},'freq_index')==1);
        freq_index = varargin{n+1};
    elseif (strcmp(varargin{n},'xaxis')==1);
        xaxis = varargin{n+1};
    elseif (strcmp(varargin{n},'param')==1);
        param = varargin{n+1};
    else
        warning('openEMS:plotFFcocx',['unknown argument key: ''' varargin{n} '''']);
    end
end

% Extract E_theta and E_phi from the nf2ff structure
E_theta = nf2ff.E_theta{freq_index} / max(nf2ff.E_norm{freq_index}(:));
E_phi   = nf2ff.E_phi{freq_index} / max(nf2ff.E_norm{freq_index}(:));

% Convert from E_theta and E_phi to co- and cross-pol as per
% Ludwig's third definition.
% Ref.: Milligan, Modern Antenna Design, 2005, p.22
% [ E_co] = [cos(phi)  -sin(phi)] [E_phi  ]
% [ E_cx] = [sin(phi)   cos(phi)] [E_theta]

E_co = E_theta.*cos(nf2ff.phi) - E_phi.*sin(nf2ff.phi);
E_cx = E_theta.*sin(nf2ff.phi) + E_phi.*cos(nf2ff.phi);

E_co_log = 20*log10(E_co) + 10*log10(nf2ff.Dmax(freq_index));
E_cx_log = 20*log10(E_cx) + 10*log10(nf2ff.Dmax(freq_index));

if (strcmp(xaxis,'theta')==1);
    xax = nf2ff.theta;
    yax1 = E_co_log(:,param);
    yax2 = E_cx_log(:,param);
    parval = nf2ff.phi(param);
    param = 'phi';
elseif (strcmp(xaxis,'phi')==1);
    xax = nf2ff.phi;
    yax = D_log(param,:);
    parval = nf2ff.theta(param);
    param = 'theta';
else
    error('openEMS:plotFFcocx','unknown parameter to ''xaxis''');
end

figure;
h = plot(xax/pi*180, yax1, xax/pi*180, yax2, '--' );

axis([-180 180]);
xlabel( sprintf('%s (deg)',xaxis ));
ylabel( 'directivity (dBi)');

createlegend = @(d)sprintf('%s = %3.1f',param,d / pi * 180);
legendtext = arrayfun(createlegend,[parval, parval],'UniformOutput',0);
legend( legendtext );
title( sprintf('far field pattern @ f = %e Hz',nf2ff.freq(freq_index)) );
grid on;

if (nargout == 0)
  clear h;
end

end
