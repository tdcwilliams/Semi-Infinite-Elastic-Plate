format long;

% paths from this repo
addpath main_funs;
addpath RTS_progs;
addpath ND_progs;
addpath compressive;

% paths from MatlabLibrary
libdir = getenv('MATLAB_LIB_DIR');
addpath([libdir '/GEN_progs']);
addpath([libdir '/SF_progs']);
addpath([libdir '/OP_progs']);
