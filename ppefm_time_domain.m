function [EEF, IEF_Ey, IEF_Ez, lt_response, time_utc] = ppefm_time_domain(Vsw, IMF_By, IMF_Bz, varargin)
% PPEFM_TIME_DOMAIN - Convert solar wind data to Equatorial Ionospheric Electric Field
%
% This function implements the Prompt Penetration Electric Field Model (PPEFM)
% using time-domain IIR filtering. It converts solar wind parameters (speed
% and IMF components) to the equatorial ionospheric eastward electric field.
%
%
% REFERENCE:
% Manoj, C., S. Maus, H. Luhr, and P. Alken (2008),
% "Comprehensive study of the dusk sector undershielding penetration electric
% field variation during high-speed solar wind streams,"
% Journal of Geophysical Research, 113, A00A17, doi:10.1029/2008JA013323
% Manoj, C., and S. Maus (2012), A real-time forecast service for the ionospheric 
% equatorial zonal electric field, Space Weather, 10, S09002, doi:10.1029/2012SW000825. 
% ==========================================================================
% SYNTAX:
%   [EEF, IEF_Ey, IEF_Ez] = ppefm_time_domain(Vsw, IMF_By, IMF_Bz)
%   [EEF, IEF_Ey, IEF_Ez, lt_response, time_utc] = ppefm_time_domain(..., 'Name', Value)
%
% ==========================================================================
% INPUTS:
%   Vsw      : Solar wind bulk speed [km/s]
%   IMF_By   : IMF Y-component in GSM [nT]
%   IMF_Bz   : IMF Z-component in GSM [nT]
%
% OPTIONAL NAME-VALUE ARGUMENTS:
%   'cadence_minutes' : Cadence in minutes (default: 5)
%   'longitude_deg'   : Geographic longitude [-180..180] or [0..360]
%                       Required if 'apply_lt' is true.
%   'start_time_utc'  : Start time of the first sample (UTC). Accepts:
%                       - datetime (timezone ignored; treated as UTC)
%                       - numeric Unix seconds since 1970-01-01
%                       Required if 'apply_lt' is true.
%   'apply_lt'        : Apply LT response modulation (default: true)
%   'apply_delay'     : Apply propagation delay to LT start time (default: true)
%   'delay_minutes'   : Propagation delay from bow shock to ionosphere (default: 17)
%   'verbose'         : Display processing info (default: false)
%
% ==========================================================================
% OUTPUTS:
%   EEF         : Equatorial ionospheric eastward electric field [mV/m]
%   IEF_Ey      : Interplanetary electric field Ey [mV/m]
%   IEF_Ez      : Interplanetary electric field Ez [mV/m]
%   lt_response : Local time response factor (length N). Empty if apply_lt=false
%   time_utc    : UTC time vector (datetime) aligned to input samples. If
%                 apply_delay=true, this is shifted forward by delay_minutes.
%
% ==========================================================================
% PHYSICAL BACKGROUND:
%   Interplanetary Electric Field (IEF): E = -V x B
%   Using GSM coordinates (solar wind ~ -X):
%     IEF_Ey = -Vsw * Bz / 1000  [mV/m]
%     IEF_Ez = -Vsw * By / 1000  [mV/m]
%
%   Unit conversion:
%     1 km/s * 1 nT = 10^3 m/s * 10^-9 T = 10^-6 V/m = 1 microV/m
%     => divide by 1000 to get mV/m
%
%   LT response is applied as:
%     EEF = LT_response * (EEF_from_Ey + EEF_from_Ez)
%
% ==========================================================================
% NOTE ON PROPAGATION DELAY:
%   OMNI data are already propagated to the bow shock nose. The software
%   applies an additional fixed 17-minute delay (Manoj et al., 2008) to account
%   for propagation from bow shock to equatorial ionosphere. This function
%   implements that same delay by shifting the LT response start time and
%   returned time vector forward by delay_minutes.
%
% ==========================================================================
% VERSION:
%   v3.1 - February 2026 - Added LT response and propagation delay matching
%          EEFMTools.GetEEF + ppefm.java logic.
%
% ==========================================================================    

    %% Parse inputs
    p = inputParser;
    addRequired(p, 'Vsw', @(x) isnumeric(x) && isvector(x));
    addRequired(p, 'IMF_By', @(x) isnumeric(x) && isvector(x));
    addRequired(p, 'IMF_Bz', @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'cadence_minutes', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'longitude_deg', [], @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'start_time_utc', [], @(x) (isnumeric(x) && isscalar(x)) || isdatetime(x));
    addParameter(p, 'apply_lt', true, @islogical);
    addParameter(p, 'apply_delay', true, @islogical);
    addParameter(p, 'delay_minutes', 17, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'verbose', false, @islogical);
    parse(p, Vsw, IMF_By, IMF_Bz, varargin{:});

    cadence_minutes = p.Results.cadence_minutes;
    longitude_deg = p.Results.longitude_deg;
    start_time_utc = p.Results.start_time_utc;
    apply_lt = p.Results.apply_lt;
    apply_delay = p.Results.apply_delay;
    delay_minutes = p.Results.delay_minutes;
    verbose = p.Results.verbose;

    %% Input validation
    Vsw = Vsw(:);
    IMF_By = IMF_By(:);
    IMF_Bz = IMF_Bz(:);

    ndata = length(Vsw);
    if length(IMF_By) ~= ndata || length(IMF_Bz) ~= ndata
        error('ppefm_time_domain:InputSizeMismatch', ...
              'All input vectors (Vsw, IMF_By, IMF_Bz) must have the same length.');
    end

    if cadence_minutes ~= 5
        warning('ppefm_time_domain:NonStandardCadence', ...
                ['Transfer function was derived for 5-minute cadence. ' ...
                 'Using %.1f-minute cadence may affect accuracy.'], cadence_minutes);
    end

    %% Step 1: Interplanetary Electric Field (IEF)
    IEF_Ey = -Vsw .* IMF_Bz / 1000;  % [mV/m]
    IEF_Ez = -Vsw .* IMF_By / 1000;  % [mV/m]

    %% Step 2: Transfer function coefficients (Dec 2008)
    Tx_a = [1.0000, -0.6023, -0.6600, 0.4451, -0.0463];
    Tx_b = [0.0052, 0.0151, 0.0014, -0.0152, -0.0061];

    Ty_a = [1.0000, -1.256299, 0.412882];
    Ty_b = [0.000824, 0.000170, -0.000296, -0.001166];

    %% Step 3: Apply IIR filter (Java CalculateTFResponse equivalent)
    EEF_from_Ey = iir_filter_time_domain(Tx_b, Tx_a, IEF_Ey);
    EEF_from_Ez = iir_filter_time_domain(Ty_b, Ty_a, IEF_Ez);

    %% Step 4: Local-time response
    lt_response = [];
    time_utc = [];
    if apply_lt
        if isempty(longitude_deg) || isempty(start_time_utc)
            error('ppefm_time_domain:LTRequiresTimeAndLongitude', ...
                  'When apply_lt is true, provide ''longitude_deg'' and ''start_time_utc''.');
        end

        % Convert start time to Unix seconds
        start_unix = to_unix_seconds(start_time_utc);

        % Apply propagation delay (bow shock -> ionosphere)
        if apply_delay && delay_minutes > 0
            start_unix = start_unix + delay_minutes * 60;
        end

        deltaT_seconds = cadence_minutes * 60;
        lt_response = calculate_lt_response(start_unix, ndata, deltaT_seconds, longitude_deg);

        % Time vector aligned with input samples, shifted by delay if applied
        time_utc = unix_to_datetime(start_unix, ndata, deltaT_seconds);
    end

    %% Step 5: Combine outputs
    if apply_lt
        EEF = lt_response(:) .* (EEF_from_Ey + EEF_from_Ez);
    else
        EEF = EEF_from_Ey + EEF_from_Ez;
    end

    if verbose
        fprintf('PPEFM Time-Domain Filter\n');
        fprintf('========================\n');
        fprintf('Input data length: %d samples\n', ndata);
        fprintf('Cadence: %.1f minutes\n', cadence_minutes);
        fprintf('IEF_Ey range: %.4f to %.4f mV/m\n', min(IEF_Ey), max(IEF_Ey));
        fprintf('IEF_Ez range: %.4f to %.4f mV/m\n', min(IEF_Ez), max(IEF_Ez));
        if apply_lt
            fprintf('LT response range: %.3f to %.3f\n', min(lt_response), max(lt_response));
            if apply_delay
                fprintf('Applied propagation delay: %.1f minutes\n', delay_minutes);
            end
        end
        fprintf('EEF range: %.4f to %.4f mV/m\n', min(EEF), max(EEF));
    end
end

%% ========================================================================
% Helper:  IIR filter 
% ========================================================================
function y = iir_filter_time_domain(b, a, x)
    x = x(:);
    b = b(:)';
    a = a(:)';

    ndata = length(x);
    nb = length(b);
    na = length(a) + 1;  

    dbuffer = zeros(nb, 1);
    y = zeros(ndata, 1);

    for j = 1:ndata
        for k = 1:(nb-1)
            dbuffer(k) = dbuffer(k+1);
        end
        dbuffer(nb) = 0.0;

        for k = 1:nb
            dbuffer(k) = dbuffer(k) + x(j) * b(k);
        end

        for k = 3:na
            dbuffer(k-1) = dbuffer(k-1) - dbuffer(1) * a(k-1);
        end

        y(j) = dbuffer(1);
    end
end

%% ========================================================================
% Helper: Local time response 
% ========================================================================
function lt_response = calculate_lt_response(start_unix, nLength, deltaT, longitude)
    % LT response lookup tables 
    lt_hours = [ ...
        -0.0833333, 0, 0.0833333, 0.1666667, 0.25, 0.3333333, ...
        0.4166667, 0.5, 0.5833333, 0.6666667, 0.75, 0.8333333, ...
        0.9166667, 1, 1.0833333, 1.1666667, 1.25, 1.3333333, ...
        1.4166667, 1.5, 1.5833333, 1.6666667, 1.75, 1.8333333, ...
        1.9166667, 2, 2.0833333, 2.1666667, 2.25, 2.3333333, ...
        2.4166667, 2.5, 2.5833333, 2.6666667, 2.75, 2.8333333, ...
        2.9166667, 3, 3.0833333, 3.1666667, 3.25, 3.3333333, ...
        3.4166667, 3.5, 3.5833333, 3.6666667, 3.75, 3.8333333, ...
        3.9166667, 4, 4.0833333, 4.1666667, 4.25, 4.3333333, ...
        4.4166667, 4.5, 4.5833333, 4.6666667, 4.75, 4.8333333, ...
        4.9166667, 5, 5.0833333, 5.1666667, 5.25, 5.3333333, ...
        5.4166667, 5.5, 5.5833333, 5.6666667, 5.75, 5.8333333, ...
        5.9166667, 6, 6.0833333, 6.1666667, 6.25, 6.3333333, ...
        6.4166667, 6.5, 6.5833333, 6.6666667, 6.75, 6.8333333, ...
        6.9166667, 7, 7.0833333, 7.1666667, 7.25, 7.3333333, ...
        7.4166667, 7.5, 7.5833333, 7.6666667, 7.75, 7.8333333, ...
        7.9166667, 8, 8.0833333, 8.1666667, 8.25, 8.3333333, ...
        8.4166667, 8.5, 8.5833333, 8.6666667, 8.75, 8.8333333, ...
        8.9166667, 9, 9.0833333, 9.1666667, 9.25, 9.3333333, ...
        9.4166667, 9.5, 9.5833333, 9.6666667, 9.75, 9.8333333, ...
        9.9166667, 10, 10.0833333, 10.1666667, 10.25, 10.3333333, ...
        10.4166667, 10.5, 10.5833333, 10.6666667, 10.75, 10.8333333, ...
        10.9166667, 11, 11.0833333, 11.1666667, 11.25, 11.3333333, ...
        11.4166667, 11.5, 11.5833333, 11.6666667, 11.75, 11.8333333, ...
        11.9166667, 12, 12.0833333, 12.1666667, 12.25, 12.3333333, ...
        12.4166667, 12.5, 12.5833333, 12.6666667, 12.75, 12.8333333, ...
        12.9166667, 13, 13.0833333, 13.1666667, 13.25, 13.3333333, ...
        13.4166667, 13.5, 13.5833333, 13.6666667, 13.75, 13.8333333, ...
        13.9166667, 14, 14.0833333, 14.1666667, 14.25, 14.3333333, ...
        14.4166667, 14.5, 14.5833333, 14.6666667, 14.75, 14.8333333, ...
        14.9166667, 15, 15.0833333, 15.1666667, 15.25, 15.3333333, ...
        15.4166667, 15.5, 15.5833333, 15.6666667, 15.75, 15.8333333, ...
        15.9166667, 16, 16.0833333, 16.1666667, 16.25, 16.3333333, ...
        16.4166667, 16.5, 16.5833333, 16.6666667, 16.75, 16.8333333, ...
        16.9166667, 17, 17.0833333, 17.1666667, 17.25, 17.3333333, ...
        17.4166667, 17.5, 17.5833333, 17.6666667, 17.75, 17.8333333, ...
        17.9166667, 18, 18.0833333, 18.1666667, 18.25, 18.3333333, ...
        18.4166667, 18.5, 18.5833333, 18.6666667, 18.75, 18.8333333, ...
        18.9166667, 19, 19.0833333, 19.1666667, 19.25, 19.3333333, ...
        19.4166667, 19.5, 19.5833333, 19.6666667, 19.75, 19.8333333, ...
        19.9166667, 20, 20.0833333, 20.1666667, 20.25, 20.3333333, ...
        20.4166667, 20.5, 20.5833333, 20.6666667, 20.75, 20.8333333, ...
        20.9166667, 21, 21.0833333, 21.1666667, 21.25, 21.3333333, ...
        21.4166667, 21.5, 21.5833333, 21.6666667, 21.75, 21.8333333, ...
        21.9166667, 22, 22.0833333, 22.1666667, 22.25, 22.3333333, ...
        22.4166667, 22.5, 22.5833333, 22.6666667, 22.75, 22.8333333, ...
        22.9166667, 23, 23.0833333, 23.1666667, 23.25, 23.3333333, ...
        23.4166667, 23.5, 23.5833333, 23.6666667, 23.75, 23.8333333, ...
        23.9166667, 24, 24.0833333 ...
    ];

    lt_values = [ ...
        -1.41069134, -1.46413114, -1.48635383, -1.50769797, -1.52815232, -1.54770561, ...
        -1.5663466, -1.58406404, -1.60084666, -1.61668322, -1.63156247, -1.64547314, ...
        -1.65840399, -1.67034376, -1.68128121, -1.69120507, -1.70010409, -1.70796703, ...
        -1.71478262, -1.72053961, -1.72523003, -1.72885896, -1.73143476, -1.73296581, ...
        -1.73346046, -1.73292708, -1.73137403, -1.72880967, -1.72524237, -1.72068049, ...
        -1.71513239, -1.70860644, -1.701111, -1.69265444, -1.68324511, -1.67289137, ...
        -1.66160161, -1.64938417, -1.63624742, -1.62219972, -1.60724944, -1.59140493, ...
        -1.57467457, -1.55706672, -1.53858974, -1.51925198, -1.49906183, -1.47802763, ...
        -1.45615776, -1.43346057, -1.40994443, -1.38561771, -1.36048875, -1.33456594, ...
        -1.30785763, -1.28037219, -1.25211907, -1.22311212, -1.19336629, -1.16289652, ...
        -1.13171774, -1.09984491, -1.06729297, -1.03407686, -1.00021153, -0.96571191, ...
        -0.93059295, -0.89486959, -0.85855678, -0.82166946, -0.78422258, -0.74623107, ...
        -0.70770987, -0.66867394, -0.62913822, -0.58911764, -0.54862716, -0.5076817, ...
        -0.46629623, -0.42448568, -0.38226499, -0.33964911, -0.29665298, -0.25329154, ...
        -0.20957974, -0.16553252, -0.12116483, -0.07649159, -0.03152777, 0.0137117, ...
        0.05921187, 0.10495781, 0.15092408, 0.19704329, 0.2432376, 0.28942912, ...
        0.33553999, 0.38149234, 0.42720831, 0.47261003, 0.51761964, 0.56215926, ...
        0.60615103, 0.64951708, 0.69217954, 0.73406056, 0.77508225, 0.81516676, ...
        0.85423622, 0.89221276, 0.92901851, 0.9645756, 0.99880618, 1.03163237, ...
        1.06297631, 1.09276012, 1.12090595, 1.14733592, 1.17197217, 1.19473684, ...
        1.21555204, 1.23433993, 1.25102262, 1.26552226, 1.27776097, 1.28766089, ...
        1.29514416, 1.3001329, 1.30258411, 1.30259422, 1.30029452, 1.2958163, ...
        1.28929086, 1.28084947, 1.27062343, 1.25874404, 1.24534258, 1.23055034, ...
        1.21449862, 1.1973187, 1.17914188, 1.16009944, 1.14032268, 1.11994288, ...
        1.09909134, 1.07789935, 1.0564982, 1.03501917, 1.01359356, 0.99235267, ...
        0.97142777, 0.95095016, 0.93105114, 0.91186198, 0.89351399, 0.87613845, ...
        0.85986665, 0.84482988, 0.83115944, 0.81898661, 0.80844269, 0.79965897, ...
        0.79276673, 0.78789726, 0.78514244, 0.7844364, 0.78567387, 0.78874958, ...
        0.79355825, 0.79999459, 0.80795334, 0.81732921, 0.82801692, 0.8399112, ...
        0.85290678, 0.86689837, 0.88178069, 0.89744847, 0.91379643, 0.9307193, ...
        0.94811179, 0.96586862, 0.98388453, 1.00205423, 1.02027244, 1.03843389, ...
        1.0564333, 1.0741654, 1.0915249, 1.10840652, 1.12470499, 1.14031504, ...
        1.15513138, 1.16904873, 1.18196183, 1.19376538, 1.20435412, 1.21362276, ...
        1.22146604, 1.22777866, 1.23247683, 1.23556264, 1.23705966, 1.23699146, ...
        1.23538159, 1.23225363, 1.22763115, 1.22153772, 1.21399689, 1.20503224, ...
        1.19466733, 1.18292574, 1.16983102, 1.15540675, 1.1396765, 1.12266382, ...
        1.10439229, 1.08488547, 1.06416694, 1.04226025, 1.01918898, 0.99497669, ...
        0.96964695, 0.94322333, 0.91572939, 0.8871887, 0.85762483, 0.82706135, ...
        0.79552182, 0.7630298, 0.72960888, 0.69528261, 0.66007455, 0.62400829, ...
        0.58710738, 0.5493954, 0.51090015, 0.47166645, 0.43174337, 0.39117996, ...
        0.35002528, 0.30832841, 0.26613839, 0.2235043, 0.18047519, 0.13710013, ...
        0.09342817, 0.04950838, 0.00538982, -0.03887844, -0.08324735, -0.12766785, ...
        -0.17209086, -0.21646734, -0.26074822, -0.30488443, -0.34882692, -0.39252662, ...
        -0.43593447, -0.47900141, -0.52167838, -0.56391631, -0.60566615, -0.64687883, ...
        -0.68750529, -0.72749648, -0.76680331, -0.80537675, -0.84316771, -0.88012715, ...
        -0.916206, -0.95135519, -0.98553617, -1.01875237, -1.0510177, -1.0823461, ...
        -1.1127515, -1.14224782, -1.17084899, -1.19856893, -1.22542158, -1.25142086, ...
        -1.27658069, -1.30091501, -1.32443775, -1.34716282, -1.36910416, -1.39027569, ...
        -1.41069134, -1.46413114, -1.48635383 ...
    ];

    % Longitude normalization 
    if longitude < -180 || longitude > 360
        longitude = 0;
    elseif longitude > 180 && longitude <= 360
        longitude = -(360 - longitude);
    end

    % UT hour from Unix seconds
    start_ut = mod(start_unix, 24*3600) / 3600;
    start_lt = start_ut + longitude / 15;
    if start_lt < 0.0
        start_lt = start_lt + 24.0;
    end

    lt_hour_array = zeros(nLength, 1);
    for i = 1:nLength
        if start_lt > 24.0
            start_lt = start_lt - 24.0;
        end
        lt_hour_array(i) = start_lt;
        start_lt = start_lt + (deltaT / 3600);
    end

    lt_response = interp1(lt_hours, lt_values, lt_hour_array, 'linear');
end

%% ========================================================================
% Helper: datetime conversions
% ========================================================================
function unix_sec = to_unix_seconds(t)
    if isdatetime(t)
        % Treat as UTC
        unix_sec = posixtime(t);
    else
        unix_sec = double(t);
    end
end

function dt = unix_to_datetime(start_unix, nLength, deltaT)
    % Returns a UTC datetime vector aligned to input samples
    dt = datetime(start_unix, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    dt = dt + seconds((0:(nLength-1)) * deltaT);
end
