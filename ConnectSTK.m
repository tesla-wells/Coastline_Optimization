%Reference 1: http://help.agi.com/stk/index.htm#training/StartMatlab.htm


clear;
%Then start up the program itself

app = actxserver('STK11.application')
app.UserControl = 1
root = app.Personality2

%Setup Scenario
scenario = root.Children.New('eScenario','MATLAB_Starter')
scenario.SetTimePeriod('Today','+24hr') %Set Animation time to like... 8 months?
root.ExecuteCommand('Animate * Reset')


%Default script for adding satellites
satellite = scenario.Children.New('eSatellite', 'LeoSat')

% IAgSatellite satellite: Satellite object
keplerian = satellite.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical'); % Use the Classical Element interface
keplerian.SizeShapeType = 'eSizeShapeAltitude';  % Changes from Ecc/Inc to Perigee/Apogee Altitude
keplerian.LocationType = 'eLocationTrueAnomaly'; % Makes sure True Anomaly is being used
keplerian.Orientation.AscNodeType = 'eAscNodeLAN'; % Use LAN instead of RAAN for data entry

% Assign the perigee and apogee altitude values:
keplerian.SizeShape.PerigeeAltitude = 500;      % km
keplerian.SizeShape.ApogeeAltitude = 600;       % km

% Assign the other desired orbital parameters:
keplerian.Orientation.Inclination = 90;         % deg
keplerian.Orientation.ArgOfPerigee = 12;        % deg
keplerian.Orientation.AscNode.Value = 24;       % deg
keplerian.Location.Value = 180;                 % deg

% Apply the changes made to the satellite's state and propagate:
satellite.Propagator.InitialState.Representation.Assign(keplerian);
satellite.Propagator.Propagate;