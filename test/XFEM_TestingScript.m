%% XFEMM Conversion
path_Lam = 'C:\GIT\xfemm\test\Lamination_1a1c094ceeed84b2.mat';
Motor = LoadDB(path_Lam);
%% write file 
fname = fullfile(tempdir, sprintf('MotorFEM_%d%d.fem', randi(1e9), 0));
fid = fopen(fname, 'w+');
fwrite(fid, Motor.FEMFile);
fclose(fid);
%% read into structure
FemmProblem = loadfemmfile(fname);
%% setup conditions
FemmProblem = setcircuitcurrent (FemmProblem, 'PhaseA', Current*exp(1i*CurAng));
FemmProblem = setcircuitcurrent (FemmProblem, 'PhaseB', Current*exp(-1i*2*pi/3+1i*CurAng));
FemmProblem = setcircuitcurrent (FemmProblem, 'PhaseC', Current*exp( 1i*2*pi/3+1i*CurAng));
FemmProblem = 

mi_modifycircprop('PhaseA', 1, Current*exp(1i*CurAng));
mi_modifycircprop('PhaseB', 1, Current*exp(-1i*2*pi/3+1i*CurAng));
mi_modifycircprop('PhaseC', 1, Current*exp( 1i*2*pi/3+1i*CurAng));
mi_modifyboundprop('periAir1', 10, p.Results.MechanicalAngle+refAngMech);
mi_modifyboundprop('periAir2', 10, -(p.Results.MechanicalAngle+refAngMech));
mi_modifyboundprop('antiAir1', 10, p.Results.MechanicalAngle+refAngMech);
mi_modifyboundprop('antiAir2', 10, -(p.Results.MechanicalAngle+refAngMech));



%%
filename = fmesher(fname);