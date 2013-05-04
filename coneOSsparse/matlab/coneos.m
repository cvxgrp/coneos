function [ varargout ] = coneos( varargin )
% coneOS 1.0
data = varargin{1};
K = varargin{2};
pars = varargin{3};
if (isfield(pars,'USE_INDIRECT') && pars.USE_INDIRECT)
    [  varargout{1:nargout}  ] = coneos_indirect( data, K, pars);
else
    [  varargout{1:nargout}  ] = coneos_direct( data, K, pars);
end