opts = { '-O', '-g' };
opts = { opts{:}, '-v' '-lm'} ;

switch computer
  case {'PCWIN', 'PCWIN64'}
    opts = {opts{:}, '-DWINDOWS'};
    opts = {opts{:}, '-ID:\libs\boost_1_57_0\'};;

  case {'GLNX86', 'GLNXA64'}
    opts = {opts{:}, '-DLINUX' };
%     opts = {opts{:}, '-I/home/mgong/libs/boost/include'};

  otherwise
    error(['Unsupported architecture ', computer, '. Please edit this M-file to fix the issue.']);
end

mex('funRegC.c', opts{:});
mex('proximalRegC.c', opts{:});



