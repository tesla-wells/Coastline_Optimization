%%% 16.346 Final Project: Coastline Optimization 
%%% Andrew Adams and Tesla Wells

%%% This code first generates points for coastlines using the mapping
%%% toolbox in MATLAB. 

%%% Then, it opens an instance of STK, given the number of coastal points
%%% to analyze, satellite parameters. Runs 180 day mission from
%%% the current time to +4320 hr.

%% Clear previous
clear all
clc
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User inputs begin here
num_points = 400;   % number of points to evaluate
% mission_start = 1;  %replace with datetime objects
% mission_end = 1;    %replace with datetime objects

% 
% % Input Orbit Parameters here for custom orbit (Comment out if using an
% % orbit below
semi = 6781000;         % semimajor axis, meters
ecc = 0.0;       % eccentricity (magnitude)
inc = 90;            % inclination, degrees
peri = 0.0;         % Argument of perigee
RAAN = 0.0;         % RAAN
true = 0.0;         % Initial true anomaly

%%%%%% Uncomment ONE subset of orbit parameters below for a specific orbit
% % Orbit parameters: ISS
% semi = 6781000;         % semimajor axis, meters
% ecc = 0.0245497;       % eccentricity (magnitude)
% inc = 51.6;            % inclination, degrees
% peri = 314.191;         % Argument of perigee
% RAAN = 306.615;         % RAAN
% true = 99.8877;         % Initial true anomaly

% % Orbit parameters, Glory
% semi = 7085000;         % semimajor axis, meters
% ecc = 0.0;       % eccentricity (magnitude)
% inc = 98.5;            % inclination, degrees
% peri = 0.0;         % Argument of perigee
% RAAN = 13.5/24*2*3.1415926;         % RAAN
% true = 0;         % Initial true anomaly
% % 
% % Orbit parameters, NPP
% semi = 7204300;         % semimajor axis, meters
% ecc = (833.7+6380)/7204.3-1;       % eccentricity (magnitude)
% inc = 98.7;            % inclination, degrees
% peri = 0.0;         % Argument of perigee
% RAAN = 0.0;         % RAAN
% true = 0.0;         % Initial true anomaly

% % Orbit parameters, USA 238
% semi = 7484000;         % semimajor axis, meters
% ecc = 1-(1049.9+6380)/7484;       % eccentricity (magnitude)
% inc = 63.4414;            % inclination, degrees
% peri = 0.0;         % Argument of perigee
% RAAN = 0.0;         % RAAN
% true = 0.0;         % Initial true anomaly

%%%%% Satellite Parameters 

% Physical satellite properties 
Sensor_FOV = 16;        % degrees

% Drag properties (if desired)
Drag = 0;               % determines whether drag model or standard propogator is used
Cd_sat =2.0;            % unitless
area_sat =0.026;        % m^2
mass_sat = 4;           % kg

%%%%%%%%%%%%%%%%%%%%%%%%%%% User inputs end here
%% Generate coastal points (Mapping toolbox, MATLAB)
load coastlines

j = 1;
ind = [0];
for i = 1:length(coastlat)
    if isnan(coastlat(i))
        k = 1;
    else
    latcopy(j,1) = coastlat(i);
    loncopy(j,1) = coastlon(i);
    j = j+1 ;
    end 
end

z = zeros(length(latcopy),1);
coastpoints = [latcopy, loncopy, z];
crash = 0;

%% STK Starter

app = actxserver('STK11.application')
app.UserControl = 1
root = app.Personality2

%% Scenario properties

scenario = root.Children.New('eScenario','ISS_J4')
scenario.SetTimePeriod('Today','+4320hr')
root.ExecuteCommand('Animate * Reset')

%% Set up target points pseudo-evenly around coastlines
n = 'CoastPoint';
points = 9600/num_points;

for i = 1:num_points
name = [n,num2str(i)];
target = scenario.Children.New('eTarget',name);
target.Position.AssignGeodetic(coastpoints(points*i,1),coastpoints(points*i,2),0)
targetarray(i)= target;
clear target
end

%% Set up satellite
satellite = scenario.Children.New('eSatellite','CubeSat');
sensor = satellite.Children.New('esensor','IRCam');
sensor.CommonTasks.SetPatternSimpleConic(Sensor_FOV,3);
satellite.MassProperties.Mass = mass_sat;

% Inertias from BeaverCube
satellite.MassProperties.Inertia.Ixx = 0.0333333;
satellite.MassProperties.Inertia.Iyy = 0.00666667;
satellite.MassProperties.Inertia.Izz = 0.0333333;

am = area_sat/mass_sat;
%% Propogate the orbit
% Set up drag model if you so desire

if Drag
% IAgSatellite satellite: Satellite object
disp(satellite.Propagator)
satellite.SetPropagatorType('ePropagatorHPOP');
set(satellite.Propagator,'Step',60);
%%
disp(satellite.Propagator.ForceModel)
endforceModel = satellite.Propagator.ForceModel;
forceModel.CentralBodyGravity.File = 'C:\Program Files\AGI\STK 11\STKData\CentralBodies\Earth\WGS84_EGM96.grv';
forceModel.CentralBodyGravity.MaxDegree = 21;
forceModel.CentralBodyGravity.MaxOrder = 21;
forceModel.Drag.Use =1;
forceModel.Drag.DragModel.Cd =2.0;
forceModel.Drag.DragModel.AreaMassRatio =0.0065;
forceModel.SolarRadiationPressure.Use =0;

integrator = satellite.Propagator.Integrator;
integrator.DoNotPropagateBelowAlt=-1e6;
integrator.IntegrationModel=3;
integrator.StepSizeControl.Method=1;
integrator.StepSizeControl.ErrorTolerance=1e-13;
integrator.StepSizeControl.MinStepSize=0.1;
integrator.StepSizeControl.MaxStepSize=30;
integrator.Interpolation.Method=1;
integrator.Interpolation.Order=7;

%%
try
root.ExecuteCommand(['SetState */Satellite/CubeSat Classical HPOP "',scenario.StartTime,'" "',scenario.StopTime,'" 60 ICRF "',scenario.StartTime,'" ',num2str(semi),' ',num2str(ecc),' ',num2str(inc),' ',num2str(peri),' ',num2str(RAAN),' ',num2str(true)])
catch ME
    switch ME.identifier
        case 'MATLAB:COM:E2147746308'
            crash = 1;
            warning('Satellite crashed into the Earth. You are not staying in space today')
        otherwise rethrow(ME)
    end

end
% or don't
else
   root.ExecuteCommand(['SetState */Satellite/CubeSat Classical J4Perturbation "',scenario.StartTime,'" "',scenario.StopTime,'" 60 ICRF "',scenario.StartTime,'" ',num2str(semi),' ',num2str(ecc),' ',num2str(inc),' ',num2str(peri),' ',num2str(RAAN),' ',num2str(true)])

end
%% Access data from all points, gather start and stop times

StartTimes = zeros(500,400);
StopTimes = zeros(500,400);
k = 1;
scenariostart = datetime(scenario.StartTime,'InputFormat','dd MMM yyyy HH:mm:ss.sss');
scenarioend = datetime(scenario.StopTime,'InputFormat','dd MMM yyyy HH:mm:ss.sss');
for i = 1:num_points
prog = i/num_points
accessarray(i) = sensor.GetAccessToObject(targetarray(i));
accessarray(i).ComputeAccess;
DParray(i) = accessarray(i).DataProviders.Item('Access Data').Exec(scenario.StartTime,scenario.StopTime);

try 
Start = DParray(i).DataSets.GetDataSetByName('Start Time').GetValues;
Stop = DParray(i).DataSets.GetDataSetByName('Stop Time').GetValues;

for j = 1:length(Start(:,1))
    Startchar = char(Start{j,1});
    Stopchar = char(Stop{j,1});
if length(Startchar) == 29 && length(Stopchar) == 29
    StartTimes(i,j) = seconds(datetime(Startchar(1:19),'InputFormat','d MMM yyyy HH:mm:ss')-scenariostart);
    StopTimes(i,j) = seconds(datetime(Stopchar(1:19),'InputFormat','d MMM yyyy HH:mm:ss')-scenariostart);
elseif length(Startchar) == 30 && length(Stopchar) == 30
    StartTimes(i,j) = seconds(datetime(Startchar(1:20),'InputFormat','dd MMM yyyy HH:mm:ss')-scenariostart);
    StopTimes(i,j) = seconds(datetime(Stopchar(1:20),'InputFormat','dd MMM yyyy HH:mm:ss')-scenariostart);
elseif length(Startchar) == 30 && length(Stopchar) == 29
    StartTimes(i,j) = seconds(datetime(Startchar(1:20),'InputFormat','dd MMM yyyy HH:mm:ss')-scenariostart);
    StopTimes(i,j) = seconds(datetime(Stopchar(1:19),'InputFormat','d MMM yyyy HH:mm:ss')-scenariostart);
elseif length(Startchar) == 29 && length(Stopchar) == 30
    StartTimes(i,j) = seconds(datetime(Startchar(1:19),'InputFormat','d MMM yyyy HH:mm:ss')-scenariostart);
    StopTimes(i,j) = seconds(datetime(Stopchar(1:20),'InputFormat','dd MMM yyyy HH:mm:ss')-scenariostart);
end

end

catch ME
    disp(ME.identifier)
    switch ME.identifier
        case 'MATLAB:COM:E2147942487'
            noflywarning = ['No Flyover for point ',num2str(i)];
            nofly(k) = i;
            k = k+1;
            warning(noflywarning);
            continue
        otherwise
            rethrow(ME)
    end
end

end

%% Estimate coastline coverage

TransitTimes = StopTimes - StartTimes;
totaltransit = sum(sum(TransitTimes));
totaltime = seconds(scenarioend-scenariostart);

NetTransit = totaltransit/totaltime;

if crash
    disp('You crashed, but while your satellite was in orbit...')
    disp(['Estimated Coastal flyover percentage = ',num2str(NetTransit*100),' percent'])
else
    disp('You made it! And while your satellite was in orbit...')
    disp(['Estimated Coastal flyover percentage = ',num2str(NetTransit*100),' percent'])
end
