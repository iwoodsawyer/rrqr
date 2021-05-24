%MAKE_RRQR Compilation of the RRQR mex-files
MATLAB_PATH = matlabroot;
COMPILE_OPTIONS = '';
v = ver('matlab');
matver = sscanf(v.Version, '%d.%d.%d')';
COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DMATLAB_VERSION=0x' sprintf('%02d%02d', matver(1), matver(2)) ];
MATLAB_VERSION = matver(1) + matver(2)/100;

if strcmpi('GLNX86', computer) || strcmpi('GLNXA64', computer) ...
        || strcmpi('MACI', computer) || strcmpi('MAC', computer) ...
        || strcmpi('MACI64', computer)
    % GNU/Linux (x86-32 or x86-64) or MacOS (Intel or PPC)
    LAPACK_PATH = ' -lmwlapack';
    if MATLAB_VERSION < 7.05
        BLAS_PATH = '';
    else
        BLAS_PATH = ' -lmwblas';
    end
    if strcmpi('GLNX86', computer) || strcmpi('GLNXA64', computer)
        COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DSkip_f2c_Undefs',' -DNON_UNIX_STDIO'];
    else
        COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DSkip_f2c_Undefs'];
    end
elseif strcmpi('PCWIN', computer) || strcmpi('PCWIN64', computer)
    % Windows (x86-32 or x86-64)
    if strcmpi('PCWIN', computer)
        if MATLAB_VERSION < 7.06
            MANUFACTURER = 'lcc';
        else
            cc = mex.getCompilerConfigurations('C','Selected');
            MANUFACTURER = cc.Manufacturer;
        end
        switch lower(MANUFACTURER)
            case {'lcc'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'lcc', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'lcc', 'libmwlapack.lib');
                COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DLCCWIN32'];
            case {'microsoft'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'microsoft', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'microsoft', 'libmwlapack.lib');
            case {'sybase'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'watcom', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'watcom', 'libmwlapack.lib');
                COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DWATCOMWIN32'];
            otherwise
                disp('Try "mex -setup", because BLAS/LAPACK library is not available!')
        end
    else
        cc = mex.getCompilerConfigurations('C','Selected');
        MANUFACTURER = cc.Manufacturer;
        switch lower(MANUFACTURER)
            case {'gnu'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'mingw64', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'mingw64', 'libmwlapack.lib');
            case {'microsoft'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
            otherwise
                disp('Try "mex -setup", because BLAS/LAPACK library is not available!')
        end
    end
    if MATLAB_VERSION < 7.05
        BLAS_PATH = ''; % On <= 7.4, BLAS in included in LAPACK
    else
        BLAS_PATH = [' "' BLAS_PATH '"'];
    end
    LAPACK_PATH = [' "' LAPACK_PATH '"'];
    RRQR_LIB = ' librrqr.lib';
    %RRQR_LIB = ' librrqr.lib libf77.lib libi77.lib'; %if you have used f2c%
    COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DMSDOS',' -DUSE_CLOCK',' -DNO_ONEXIT'];
else
    error('Unsupported platform')
end

% Large array dims for 64 bits platforms appeared in Matlab 7.3
if (strcmpi('GLNXA64', computer) || strcmpi('PCWIN64', computer) ...
        || strcmpi('MACI64', computer))
    if ~(MATLAB_VERSION < 9.04)
        COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -R2018a' ];
    elseif ~(MATLAB_VERSION < 7.03)
        COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -largeArrayDims' ];
    end
end

% Comment next line to suppress optimization
COMPILE_OPTIONS = [ ' -O' COMPILE_OPTIONS ];

% Comment next line to suppress compilation debugging info
%COMPILE_OPTIONS = [ ' -v' COMPILE_OPTIONS ];

disp('Compiling rrqrx...')
eval(['mex ', COMPILE_OPTIONS, ' -I.', ' -Ilibf2c', ' rrqrx.c',...
    ' acm782/zgeqpb.c',' acm782/dgeqpb.c',' acm782/cgeqpb.c',' acm782/sgeqpb.c',...
    ' acm782/zgeqpc.c',' acm782/dgeqpc.c',' acm782/cgeqpc.c',' acm782/sgeqpc.c',...
    ' acm782/zgeqpw.c',' acm782/dgeqpw.c',' acm782/cgeqpw.c',' acm782/sgeqpw.c',...
    ' acm782/zgeqpx.c',' acm782/dgeqpx.c',' acm782/cgeqpx.c',' acm782/sgeqpx.c',...
    ' acm782/ztrqpx.c',' acm782/dtrqpx.c',' acm782/ctrqpx.c',' acm782/strqpx.c',...
    ' acm782/ztrqxc.c',' acm782/dtrqxc.c',' acm782/ctrqxc.c',' acm782/strqxc.c',...
    ' acm782/zlasmx.c',' acm782/dlasmx.c',' acm782/clasmx.c',' acm782/slasmx.c',...
    ' acm782/zlauc1.c',' acm782/dlauc1.c',' acm782/clauc1.c',' acm782/slauc1.c',...
    ' acm782/ztrrnk.c',' acm782/dtrrnk.c',' acm782/ctrrnk.c',' acm782/strrnk.c',...
    ' libf2c/d_cnjg.c',' libf2c/r_cnjg.c',' libf2c/d_sign.c',' libf2c/r_sign.c',...
    ' libf2c/z_abs.c',' libf2c/c_abs.c',' libf2c/cabs.c',' libf2c/pow_dd.c',' libf2c/s_copy.c ',BLAS_PATH, LAPACK_PATH]);
disp('Compiling rrqry...')
eval(['mex ', COMPILE_OPTIONS, ' -I.', ' -Ilibf2c', ' rrqry.c',...
    ' acm782/zgeqpb.c',' acm782/dgeqpb.c',' acm782/cgeqpb.c',' acm782/sgeqpb.c',...
    ' acm782/zgeqpc.c',' acm782/dgeqpc.c',' acm782/cgeqpc.c',' acm782/sgeqpc.c',...
    ' acm782/zgeqpw.c',' acm782/dgeqpw.c',' acm782/cgeqpw.c',' acm782/sgeqpw.c',...
    ' acm782/zgeqpy.c',' acm782/dgeqpy.c',' acm782/cgeqpy.c',' acm782/sgeqpy.c',...
    ' acm782/ztrqpy.c',' acm782/dtrqpy.c',' acm782/ctrqpy.c',' acm782/strqpy.c',...
    ' acm782/ztrqyc.c',' acm782/dtrqyc.c',' acm782/ctrqyc.c',' acm782/strqyc.c',...
    ' acm782/ztrqxc.c',' acm782/dtrqxc.c',' acm782/ctrqxc.c',' acm782/strqxc.c',...
    ' acm782/zlasmx.c',' acm782/dlasmx.c',' acm782/clasmx.c',' acm782/slasmx.c',...
    ' acm782/zlauc1.c',' acm782/dlauc1.c',' acm782/clauc1.c',' acm782/slauc1.c',...
    ' acm782/ztrrnk.c',' acm782/dtrrnk.c',' acm782/ctrrnk.c',' acm782/strrnk.c',...
    ' libf2c/d_cnjg.c',' libf2c/r_cnjg.c',' libf2c/d_sign.c',' libf2c/r_sign.c',...
    ' libf2c/z_abs.c',' libf2c/c_abs.c',' libf2c/cabs.c',' libf2c/pow_dd.c',' libf2c/s_copy.c ',BLAS_PATH, LAPACK_PATH]);
disp('Compilation of mex files was successful.')