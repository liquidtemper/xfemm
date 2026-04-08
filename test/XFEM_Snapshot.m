function [snapshot, mesh, FEMFile, ansFile] = XFEM_Snapshot(Motor, varargin)
p=inputParser;
p.addParameter('MechanicalAngle', 0);
p.addParameter('CurrentAngle', 0);
p.addParameter('CurrentDensity', 0);
p.addParameter('ElectricFrequency', 0);
p.addParameter('MagnetTemperature', 20);
p.addParameter('FEMNum', randi(1e9));
p.addParameter('keepOpen', false);
p.addParameter('getMesh', false);
p.addParameter('prevSol', '');
p.parse(varargin{:});

t = getCurrentTask(); 
if isempty(t); t.ID = 1;end
tempdir = fullfile('C:', 'femm42', 'Temp');
if ~exist(tempdir, 'file'); mkdir(tempdir); end
fname = fullfile(tempdir, sprintf('MotorFEM_%d%d.fem', p.Results.FEMNum, t.ID));
if ~isfield(Motor, 'FEMFile')
    FEMFile = geom2femm(Motor.geom);
else
    fid = fopen(fname, 'w+');
    fwrite(fid, Motor.FEMFile);
    fclose(fid);
    FEMFile = Motor.FEMFile;
    openfemm(1);
    opendocument(fname);
end
mi_saveas(fname);
hideconsole;
hidepointprops;
main_minimize;

electricAngle = p.Results.MechanicalAngle * pi/180 * Motor.Configuration.Rotor.Poles/2;
% set circuit properties
CurAng = p.Results.CurrentAngle*pi/180 + electricAngle;
snapshot.Theta = electricAngle;
Current = p.Results.CurrentDensity * Motor.Configuration.Winding.SlotSurface;

if isfield(Motor, 'Parameters')
    %reference mechanical angle to center the phases. convert to mechanical
    %angles
    refAngMech = -Motor.Parameters.ReferenceAngle*2/Motor.Configuration.Rotor.Poles; 
else
    refAngMech = 0; %first snapshot
end

mi_modifycircprop('PhaseA', 1, Current*exp(1i*CurAng));
mi_modifycircprop('PhaseB', 1, Current*exp(-1i*2*pi/3+1i*CurAng));
mi_modifycircprop('PhaseC', 1, Current*exp( 1i*2*pi/3+1i*CurAng));
mi_modifyboundprop('periAir1', 10, p.Results.MechanicalAngle+refAngMech);
mi_modifyboundprop('periAir2', 10, -(p.Results.MechanicalAngle+refAngMech));
mi_modifyboundprop('antiAir1', 10, p.Results.MechanicalAngle+refAngMech);
mi_modifyboundprop('antiAir2', 10, -(p.Results.MechanicalAngle+refAngMech));
% modify rotor magnet temperature
magnetInd = find([Motor.geom.materials.aBr] ~=0);
for mm=1:length(magnetInd)
    ind     = magnetInd(mm);
    thismag = Motor.geom.materials(ind);
    T       = p.Results.MagnetTemperature;
    Hcnew   = thismag.H_c * (1+thismag.aBr/100*(T-20));
    mi_modifymaterial(thismag.name, 3,Hcnew);
end

try
    mi_minimize;
    mi_createmesh;
catch 
    pause(1)
    mi_minimize;
    mi_createmesh;
end
mi_probdef(p.Results.ElectricFrequency,'millimeters','planar',1e-8, 1, 15, 0);
if p.Results.ElectricFrequency > 0
    mi_modifymaterial('Copper', 9, 3);
end
if ~isempty(p.Results.prevSol)
    mi_setprevious(p.Results.prevSol)
end
try 
    mi_analyze;
    mi_loadsolution;
catch ME
    delete(fname);
    delete(strrep(fname, '.fem', '.ans'));
    rethrow(ME);
end

ansFilePath = regexprep(fname, '\.fem$', '.ans');

if nargout>=4
    ansFile     = fileread(ansFilePath);
end

if ~isfield(Motor, 'Mesh')
    elem = mo_getelement(1);
    elems = NaN(length(elem), mo_numelements);
    for mm=1:mo_numelements
        elems(:,mm) = mo_getelement(mm);
    end
    mesh.elements = elems;

    node = mo_getnode(1);
    nodes = NaN(length(node), mo_numnodes);
    for nn = 1:mo_numnodes
        nodes(:,nn) = mo_getnode(nn);
    end
    mesh.nodes = nodes;
else
    mesh = Motor.Mesh;
end
%% Load stator and rotor mesh elements
elem = mesh.elements(:,1);
elemprops = mo_getpointvalues(elem(4), elem(5));
blocks = mesh.elements(7,:);
[~,touse] = ismember(blocks, [21:100000]);
touse = touse>0; %only get element properties from certain block labels
snapshot.MeshValues = NaN(length(elemprops), size(mesh.elements,2));
if p.Results.getMesh
    for mm= 1:size(mesh.elements,2)
        if ~touse(mm); continue; end
        elem = mesh.elements(:,mm);
        snapshot.MeshValues(:,mm) = mo_getpointvalues(elem(4), elem(5));
    end
end
%% load rotor cage MVP from block labels
if strcmp(Motor.Configuration.MachineType, 'IM')
    snapshot.RotorCage = NaN(4,length(Motor.geom.cage.positions));
    for bb = 1:length(Motor.geom.cage.positions)
        pos = Motor.geom.cage.positions(bb);
        mo_selectblock(real(pos), imag(pos));
        snapshot.RotorCage(1,bb) = mo_blockintegral(1); %Integral of A
        snapshot.RotorCage(2,bb) = mo_blockintegral(5); %Surface Area
        snapshot.RotorCage(3,bb) = mo_blockintegral(7); %Current
        % calculate mvp on each ba
        snapshot.RotorCage(4,bb) = snapshot.RotorCage(1,bb)/snapshot.RotorCage(2,bb)*1e3;
        mo_clearblock;
    end
else
    snapshot.RotorCage = [];
end
%% get circuit properties
tmp = mo_getcircuitproperties('PhaseA');
Circuit = repmat(tmp,3,1);
Circuit(2,:) = mo_getcircuitproperties('PhaseB');
Circuit(3,:) = mo_getcircuitproperties('PhaseC');
%%
% get to single strand MVP here instead of series wound
% turns for this femm winding snapshot. Current also total current
% Neff = Motor.Configuration.Winding.Layers * Motor.Configuration.Stator.Slots/Motor.Configuration.Winding.Phases;
% snapshot.Circuit(:,3) = snapshot.Circuit(:,3)*1e3/Neff; 
phidq               = abc2dq(Circuit(:,3), electricAngle*180/pi);
curdq               = abc2dq(Circuit(:,1), electricAngle*180/pi);
snapshot.CurrentDQ0 = curdq; 
snapshot.FluxDQ0    = phidq; 
%% Get Torque
trq = mo_gapintegral('periAir1', 0);
trq = trq + mo_gapintegral('periAir2', 0);
trq = trq + mo_gapintegral('antiAir1', 0);
trq = trq + mo_gapintegral('antiAir2', 0);
snapshot.Torque = trq;
%% Airgap parameters
angle             = linspace(0,360,1081);
angle             = angle(1:end-1); %skip duplicate point
snapshot.MvpAir   = mo_getgapa('periAir1', angle);
snapshot.BAir     = mo_getgapb('periAir1', angle);
snapshot.Temperature = p.Results.MagnetTemperature;
%% Close snapshot
if any(strcmp(p.UsingDefaults, 'FEMNum')) && ~p.Results.keepOpen
    closefemm;
    delete(fname);
    delete(strrep(fname, '.fem', '.ans'));
end


end