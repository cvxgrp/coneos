function shim = cvx_coneos( shim )

% CVX_SOLVER_SHIM	SeDuMi interface for CVX.
%   This procedure returns a 'shim': a structure containing the necessary
%   information CVX needs to use this solver in its modeling framework.

if ~isempty( shim.solve ),
    return
end
if isempty( shim.name ),
    fname = 'coneos_matlab.m';
    [ fs, ps, int_path ] = cvx_version;
    int_path(end+1) = fs;
    int_plen = length( int_path );
    shim.name = 'coneos_matlab';
    shim.dualize = true;
    flen = length(fname);
    fpaths = { [ int_path, 'coneos_matlab', fs, fname ] };
    fpaths = [ fpaths ; which( fname, '-all' ) ];
    old_dir = pwd;
    oshim = shim;
    shim = [];
    for k = 1 : length(fpaths),
        fpath = fpaths{k};
        if ~exist( fpath, 'file' ) || any( strcmp( fpath, fpaths(1:k-1) ) ),
            continue
        end
        new_dir = fpath(1:end-flen-1);
        cd( new_dir );
        tshim = oshim;
        tshim.fullpath = fpath;
        tshim.version = 'unknown';
        is_internal = strncmp( new_dir, int_path, int_plen );
        if is_internal,
            tshim.location = [ '{cvx}', new_dir(int_plen:end) ];
        else
            tshim.location = new_dir;
        end
        try
            fid = fopen(fname);
            otp = fread(fid,Inf,'uint8=>char')';
            fclose(fid);
        catch errmsg
            tshim.error = sprintf( 'Unexpected error:\n%s\n', errmsg.message );
        end
        if isempty( tshim.error ),
            otp = regexp( otp, 'coneOS \d+\.\d+', 'match' );
            if ~isempty(otp), tshim.version = otp{1}(8:end); end
            tshim.check = @check;
            tshim.solve = @solve;
            if k ~= 2,
                tshim.path = [ new_dir, ps ];
            end
        end
        shim = [ shim, tshim ]; %#ok
    end
    cd( old_dir );
    if isempty( shim ),
        shim = oshim;
        shim.error = 'Could not find a coneOS installation.';
    end
else
    shim.check = @check;
    shim.solve = @solve;
end

function found_bad = check( nonls ) %#ok
found_bad = false;

function [ x, status, tol, iters, y, z ] = solve( At, b, c, nonls, quiet, prec, settings )

n = length( c );
m = length( b );
K = struct( 'f', 0, 'l', 0, 'q', [], 's', [], 'scomplex', [], 'ycomplex', [] );
reord = struct( 'n', 0, 'r', [], 'c', [], 'v', [] );
reord = struct( 'f', reord, 'l', reord, 'a', reord, 'q', reord, 's', reord, 'h', reord );
reord.f.n = n;
for k = 1 : length( nonls ),
    temp = nonls( k ).indices;
    nn = size( temp, 1 );
    nv = size( temp, 2 );
    nnv = nn * nv;
    tt = nonls( k ).type;
    reord.f.n = reord.f.n - nnv;
    if strncmp( tt, 'i_', 2 ),
        error( 'coneOS does not support integer variables.' );
    elseif nn == 1 || isequal( tt, 'nonnegative' ),
        reord.l.r = [ reord.l.r ; temp(:) ];
        reord.l.c = [ reord.l.c ; reord.l.n + ( 1 : nnv )' ];
        reord.l.v = [ reord.l.v ; ones( nnv, 1 ) ];
        reord.l.n = reord.l.n + nnv;
    elseif isequal( tt, 'lorentz' ),
        if nn == 2,
            rr = [ temp ; temp ];
            cc = reshape( floor( 1 : 0.5 : 2 * nv + 0.5 ), 4, nv );
            vv = [1;1;-1;1]; vv = vv(:,ones(1,nv));
            reord.a.r = [ reord.a.r ; rr(:) ];
            reord.a.c = [ reord.a.c ; cc(:) + reord.a.n ];
            reord.a.v = [ reord.a.v ; vv(:) ];
            reord.a.n = reord.a.n + nnv;
        else
            temp = temp( [ end, 1 : end - 1 ], : );
            reord.q.r = [ reord.q.r ; temp(:) ];
            reord.q.c = [ reord.q.c ; reord.q.n + ( 1 : nnv )' ];
            reord.q.v = [ reord.q.v ; ones(nnv,1) ];
            reord.q.n = reord.q.n + nnv;
            K.q = [ K.q, nn * ones( 1, nv ) ];
        end
    elseif isequal( tt, 'semidefinite' ),
        nn = 0.5 * ( sqrt( 8 * nn + 1 ) - 1 );
        str = cvx_create_structure( [ nn, nn, nv ], 'symmetric' );
        K.s = [ K.s, nn * ones( 1, nv ) ];
        [ cc, rr, vv ] = find( cvx_invert_structure( str, 'compact' ) );
        rr = temp( rr );
        reord.s.r = [ reord.s.r; rr( : ) ];
        reord.s.c = [ reord.s.c; cc( : ) + reord.s.n ];
        reord.s.v = [ reord.s.v; vv( : ) ];
        reord.s.n = reord.s.n + nn * nn * nv;
    elseif isequal( tt, 'hermitian-semidefinite' ) && 1,
        % SeDuMi's complex SDP support was broken with the 1.3 update. It
        % simply applies a conversion to real SDP, however, which we can
        % reproduce here.
        %   X >= 0 <==> exists [ Y1, Y2^T ; Y2, Y3 ] >= 0 s.t.
        %               Y1 + Y3 == real(X), Y2 - Y2^T == imag(X)
        nsq = nn; nn = sqrt( nn );
        str = cvx_create_structure( [ nn, nn, nv ], 'hermitian' );
        [ cc, rr, vv ] = find( cvx_invert_structure( str, 'compact' ) );
        cc = cc - 1;
        mm = floor( cc / nsq );
        cc = cc - mm * nsq;
        jj = floor( cc / nn );
        ii = cc - jj * nn + 1;
        jj = jj + 1;
        mm = mm + 1;
        vr = real( vv );
        vi = imag( vv );
        ii = [ ii + nn * ~vr ; ii + nn * ~vi ];
        jj = [ jj ; jj + nn ]; %#ok
        vv = sqrt( 0.5 ) * [ vr + vi ; vr - vi ];
        rr = [ rr ; rr ]; %#ok
        mm = [ mm ; mm ]; %#ok
        [ jj, ii ] = deal( min( ii, jj ), max( ii, jj ) );
        cc = ii + ( jj - 1 ) * ( 2 * nn ) + ( mm - 1 ) * ( 4 * nsq );
        K.s = [ K.s, 2 * nn * ones( 1, nv ) ];
        rr = temp( rr );
        reord.s.r = [ reord.s.r; rr( : ) ];
        reord.s.c = [ reord.s.c; cc( : ) + reord.s.n ];
        reord.s.v = [ reord.s.v; vv( : ) ];
        reord.s.n = reord.s.n + 4 * nsq * nv;
    elseif isequal( tt, 'hermitian-semidefinite' ),
        % If SeDuMi's complex SDP support is restored, we can usse this
        % simpler mapping instead, if we choose.
        K.scomplex = [ K.scomplex, length( K.s ) + ( 1 : nv ) ];
        nn = sqrt( nn );
        str = cvx_create_structure( [ nn, nn, nv ], 'hermitian' );
        K.s = [ K.s, nn * ones( 1, nv ) ];
        stri = cvx_invert_structure( str, 'compact' )';
        [ rr, cc, vv ] = find( stri );
        rr = temp( rr );
        reord.s.r = [ reord.s.r; rr( : ) ];
        reord.s.c = [ reord.s.c; cc( : ) + reord.s.n ];
        reord.s.v = [ reord.s.v; vv( : ) ];
        reord.s.n = reord.s.n + size( stri, 2 );
    else
        error( 'Unsupported nonlinearity: %s', tt );
    end
end
if reord.f.n > 0,
    reord.f.r = ( 1 : n )';
    reord.f.r( [ reord.l.r ; reord.a.r ; reord.q.r ; reord.s.r ] ) = [];
    reord.f.c = ( 1 : reord.f.n )';
    reord.f.v = ones(reord.f.n,1);
end
n_d = max( m - n - reord.f.n + 1, isempty( At ) );
if n_d,
    reord.l.n = reord.l.n + n_d;
end
K.f = reord.f.n;
K.l = reord.l.n + reord.a.n;
n_out = reord.f.n;
reord.l.c = reord.l.c + n_out; n_out = n_out + reord.l.n;
reord.a.c = reord.a.c + n_out; n_out = n_out + reord.a.n;
reord.q.c = reord.q.c + n_out; n_out = n_out + reord.q.n;
reord.s.c = reord.s.c + n_out; n_out = n_out + reord.s.n;
reord = sparse( ...
    [ reord.f.r ; reord.l.r ; reord.a.r ; reord.q.r ; reord.s.r ], ...
    [ reord.f.c ; reord.l.c ; reord.a.c ; reord.q.c ; reord.s.c ], ...
    [ reord.f.v ; reord.l.v ; reord.a.v ; reord.q.v ; reord.s.v ], ...
    n, n_out );

At = reord' * At;
c  = reord' * c;
pars.free = K.f > 1 && nnz( K.q );
pars.eps     = prec(1);
pars.bigeps  = prec(3);
if quiet,
    pars.fid = 0;
end
add_row = isempty( At );
if add_row,
    K.f = K.f + 1;
    At = sparse( 1, 1, 1, n_out + 1, 1 );
    b = 1;
    c = [ 0 ; c ];
end

if( ~isreal(At) || ~isreal(c) || ~isreal(b) )
    error('coneOS does not handle complex data');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.A = sparse(At);
data.b = full(c);
data.c = -full(b);

%{
% symmetrize SD cone matrices:
[mm,nn]=size(data.A);
idx = K.f + K.l + sum(K.q);
for i=1:size(K.s)
    for j=1:nn
        work = data.A(idx+1:idx+K.s(i)^2, j);
        work = reshape(work,K.s(i),K.s(i));
        if any(any(work~=work'))
            %fprintf('warning: symmetrizing A\n')
            work = (work+work')/2;
            data.A(idx+1:idx+K.s(i)^2, j) = work(:);
        end
    end
    
    work = data.b(idx+1:idx+K.s(i)^2);
    work = reshape(work,K.s(i),K.s(i));
    if any(any(work~=work'))
        %fprintf('warning: symmetrizing b\n')
        work = (work+work')/2;
        data.b(idx+1:idx+K.s(i)^2) = work(:);
    end
    
    idx = idx + K.s(i)^2;
end
clear work;
%}

[m1,n1] = size(K.q);
if m1 == 1,
    K.q = K.q';
end
[m1,n1] = size(K.s);
if m1 == 1,
    K.s = K.s';
end

if (isfield(settings,'UNDET_TOL'))
    pars.UNDET_TOL = settings.UNDET_TOL;
end
if prec(1)==0
    % just run until MAX_ITERS
    pars.EPS_ABS = 0;
end
if (isfield(settings,'EPS'))
    pars.EPS_ABS = settings.EPS;
end
if (isfield(settings,'MAX_ITERS'))
    pars.MAX_ITERS = settings.MAX_ITERS;
end
if ~quiet
    pars.VERBOSE = 1;
end
if (isfield(settings,'USE_INDIRECT'))
    pars.USE_INDIRECT = settings.USE_INDIRECT;
end
if (isfield(settings,'NORMALIZE'))
    pars.NORMALIZE = settings.NORMALIZE;
end
if (isfield(settings,'CG_TOL'))
    pars.CG_TOL = settings.CG_TOL;
end
if (isfield(settings,'RHO_X'))
    pars.RHO_X = settings.RHO_X;
end
if (isfield(settings,'CG_MAX_ITS'))
    pars.CG_MAX_ITS = settings.CG_MAX_ITS;
end
if (isfield(settings,'ALPHA'))
    pars.ALPHA = settings.ALPHA;
end
if (isfield(settings,'RHO_X'))
    pars.RHO_X = settings.RHO_X;
end


[ yy, xx, info ] = cvx_run_solver( @coneos_matlab, data, K, pars, 'xx', 'yy', 'info', settings, 3);

if add_row,
    xx = xx(2:end);
    yy = zeros(0,1);
    At = zeros(n_out,0);
    % b  = zeros(0,1);
    c  = c(2:end);
end

xx = full( xx );
yy = full( yy );
x = real( reord * xx );
y = yy;
z = real( reord * ( c - At * yy ) );
if add_row, y = zeros( 0, 1 ); end

tol = -1;%max(info.presid,info.dresid);
iters = info.iter;
status = info.status;
% coneOS targets the dual to sedumi formulation:
if (strcmp(status,'Unbounded'))
    status = 'Infeasible';
elseif (strcmp(status,'Infeasible'))
    status = 'Unbounded';
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
