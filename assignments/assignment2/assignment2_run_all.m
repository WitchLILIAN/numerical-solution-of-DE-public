% Run all MATLAB programs for Assignment 2.
clear; close all; clc;

rootDir = fileparts(mfilename('fullpath'));
set(0, 'DefaultFigureVisible', 'on');

fprintf('Running Assignment 2 Problem 1...\n');
run(fullfile(rootDir, 'assignment2_problem1.m'));

fprintf('Running Assignment 2 Problem 3...\n');
run(fullfile(rootDir, 'assignment2_problem3.m'));

fprintf('Running Assignment 2 Problem 4...\n');
run(fullfile(rootDir, 'assignment2_problem4.m'));

fprintf('\nAll computations finished. Tables were printed in the Command Window.\n');
fprintf('Figures were created with plot/mesh commands and are visible in MATLAB.\n');
