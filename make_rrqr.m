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
            cc = mex.getCompilerConfigurations('Any','Selected');
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
        cc = mex.getCompilerConfigurations('Any','Selected');
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
eval(['mex ', COMPILE_OPTIONS, ' rrqrx.c',...
    ' src/zgeqpb.c',' src/dgeqpb.c',' src/cgeqpb.c',' src/sgeqpb.c',...
    ' src/zgeqpc.c',' src/dgeqpc.c',' src/cgeqpc.c',' src/sgeqpc.c',...
    ' src/zgeqpw.c',' src/dgeqpw.c',' src/cgeqpw.c',' src/sgeqpw.c',...
    ' src/zgeqpx.c',' src/dgeqpx.c',' src/cgeqpx.c',' src/sgeqpx.c',...
    ' src/ztrqpx.c',' src/dtrqpx.c',' src/ctrqpx.c',' src/strqpx.c',...
    ' src/ztrqxc.c',' src/dtrqxc.c',' src/ctrqxc.c',' src/strqxc.c',...
    ' src/zlasmx.c',' src/dlasmx.c',' src/clasmx.c',' src/slasmx.c',...
    ' src/zlauc1.c',' src/dlauc1.c',' src/clauc1.c',' src/slauc1.c',...
    ' src/ztrrnk.c',' src/dtrrnk.c',' src/ctrrnk.c',' src/strrnk.c',...
    ' src/d_cnjg.c',' src/r_cnjg.c',' src/d_sign.c',' src/r_sign.c',...
    ' src/z_abs.c',' src/c_abs.c',' src/cabs.c',' src/pow_dd.c',' src/s_copy.c ',BLAS_PATH, LAPACK_PATH]);
disp('Compiling rrqry...')
eval(['mex ', COMPILE_OPTIONS, ' rrqry.c',...
    ' src/zgeqpb.c',' src/dgeqpb.c',' src/cgeqpb.c',' src/sgeqpb.c',...
    ' src/zgeqpc.c',' src/dgeqpc.c',' src/cgeqpc.c',' src/sgeqpc.c',...
    ' src/zgeqpw.c',' src/dgeqpw.c',' src/cgeqpw.c',' src/sgeqpw.c',...
    ' src/zgeqpy.c',' src/dgeqpy.c',' src/cgeqpy.c',' src/sgeqpy.c',...
    ' src/ztrqpy.c',' src/dtrqpy.c',' src/ctrqpy.c',' src/strqpy.c',...
    ' src/ztrqyc.c',' src/dtrqyc.c',' src/ctrqyc.c',' src/strqyc.c',...
    ' src/ztrqxc.c',' src/dtrqxc.c',' src/ctrqxc.c',' src/strqxc.c',...
    ' src/zlasmx.c',' src/dlasmx.c',' src/clasmx.c',' src/slasmx.c',...
    ' src/zlauc1.c',' src/dlauc1.c',' src/clauc1.c',' src/slauc1.c',...
    ' src/ztrrnk.c',' src/dtrrnk.c',' src/ctrrnk.c',' src/strrnk.c',...
    ' src/d_cnjg.c',' src/r_cnjg.c',' src/d_sign.c',' src/r_sign.c',...
    ' src/z_abs.c',' src/c_abs.c',' src/cabs.c',' src/pow_dd.c',' src/s_copy.c ',BLAS_PATH, LAPACK_PATH]);
disp('Compilation of mex files was successful.')