% mess_reset clears the work space variables, closes the open figures and
% clears the terminal.
%
% Note that only variable data is cleared, such that JIT compiler data can
% be preserved, which should speedup things over doing a manual 'clear all'.

close all;

if exist('OCTAVE_VERSION', 'builtin')
    clear -variables;
else
    clearvars;
end

clc;