function grad = read_siemens_grad_file(grad_file)

%% Read a Siemens .grad file
fid = fopen(grad_file);

%--------------------------------------------------------------------------
% Skip lines with # sign
%--------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    if isempty(tline) || strcmp(tline(1), '#')
        continue;
    else
        break;
    end
end

%--------------------------------------------------------------------------
% Get header information
% ' flagship 097, Gradientsystem Verio , Gx,y,z = 45 mT/m'
%--------------------------------------------------------------------------
info = tline;

equal_loc = strfind(info, '=');
slash_loc = strfind(info, '/');

if length(slash_loc) > 1
    Gx = str2double(info(equal_loc+1:slash_loc(end-1)-1));
    Gz = str2double(info(slash_loc(end-1)+1:slash_loc(end)-3));
    Gy = Gx;
else
    Gx = str2double(info(equal_loc+1:slash_loc-3));
    Gy = Gx;
    Gz = Gx;
end

%--------------------------------------------------------------------------
% win_low = 0, win_high = 0, win_algo = 0, win_dummy = 2;
%--------------------------------------------------------------------------
tline = fgetl(fid);

equal_loc = strfind(tline, '=');
comma_loc = strfind(tline, ',');

win_low   = str2double(tline(equal_loc(1)+1:comma_loc(1)-1));
win_high  = str2double(tline(equal_loc(2)+1:comma_loc(2)-1));
win_algo  = str2double(tline(equal_loc(3)+1:comma_loc(3)-1));
win_dummy = str2double(tline(equal_loc(4)+1:end-1));

%--------------------------------------------------------------------------
% 0.250 m = R0, lnorm = 4? A(1,0) = B(1,1) = A(1,1) = 0;
% Get the radious R0 of a gradient coil
%--------------------------------------------------------------------------
tline = fgetl(fid);

m_loc = strfind(tline, 'm');
equal_loc = strfind(tline, '=');
question_loc = strfind(tline, '?');

R0 = str2double(tline(1:m_loc(1)-1)); % [m]
lnorm = str2double(tline(equal_loc(2)+1:question_loc-1)); % [mT/m]

%--------------------------------------------------------------------------
% Skip 6 lines
%--------------------------------------------------------------------------
tline = fgetl(fid); % 0 = CoSyMode
tline = fgetl(fid); % Expansion in Spherical Harmonics
tline = fgetl(fid); % ================================
tline = fgetl(fid); %
tline = fgetl(fid); % NO.  TYPE             SPECTRUM                 AXIS
tline = fgetl(fid); %

%% Get Legendre coefficients in spherical harmonics
largest_order = 0;
while ~feof(fid)
    tline = fgetl(fid);
    tline = regexprep(tline, '\t', ' ');
    tline = regexprep(tline, ' +', ' ');
    if length(tline) == 1 || isempty(tline)
        break;
    end
    if isequal(tline(end), ' ') % remove an empty space at the end of a string
        tline = tline(1:end-1);
    end
    space_loc = strfind(tline, ' ');
    left_round_loc = strfind(tline, '(');
    right_round_loc = strfind(tline, ')');
    comma_loc = strfind(tline, ',');

    %----------------------------------------------------------------------
    % Parse each line
    % ' 1 A( 3, 0) -0.1073 z'
    % '101 A( 3, 1) -0.1168 x'
    %----------------------------------------------------------------------
    no       = str2double(tline(1:left_round_loc-2)); % 1 or 101 or 201
    type     = tline(left_round_loc-1); % A or B
    order    = str2double(tline(left_round_loc+1:comma_loc-1)); % ( 3, 0) or (14, 1)
    degree   = str2double(tline(comma_loc+1:right_round_loc-1));
    spectrum = str2double(tline(right_round_loc+1:space_loc(end)-1));
    axis     = tline(space_loc(end)+1);
    %fprintf('%s\n%2d %s(%2d,%2d) %6.4f %s\n', tline, no, type, order, degree, spectrum, axis);

    %pause;
    if order > largest_order
        largest_order = order;
    end

    if strcmp(type, 'A')
        switch axis
            case 'z'
                alpha_z(order+1, degree+1) = spectrum;
            case 'x'
                alpha_x(order+1, degree+1) = spectrum;
            case 'y'
                alpha_y(order+1, degree+1) = spectrum;
            otherwise
        end
    elseif strcmp(type, 'B')
        switch axis
            case 'z'
                beta_z(order+1, degree+1) = spectrum;
            case 'x'
                beta_x(order+1, degree+1) = spectrum;
            case 'y'
                beta_y(order+1, degree+1) = spectrum;
            otherwise
        end
    end
end
fclose(fid);

%% Resize the alpha and beta
if size(alpha_z,2) ~= largest_order+1
    alpha_z(largest_order+1, largest_order+1) = 0;
end

if size(alpha_x,2) ~= largest_order+1
    alpha_x(largest_order+1, largest_order+1) = 0;
end

if size(beta_y,2) ~= largest_order+1
    beta_y(largest_order+1, largest_order+1) = 0;
end

%% Define nonexisting variables as zero arrays
if ~exist('beta_z', 'var')
    beta_z = alpha_z * 0;
else
    if size(beta_z,2) ~= largest_order+1
        beta_z(largest_order+1, largest_order+1) = 0;
    end
end

if ~exist('beta_x', 'var')
    beta_x =alpha_x * 0;
else
    if size(beta_x,2) ~= largest_order+1
        beta_x(largest_order+1, largest_order+1) = 0;
    end
end
    
if ~exist('alpha_y', 'var')
    alpha_y = beta_y * 0;
else
    if size(alpha_y,2) ~= largest_order+1
        alpha_y(largest_order+1, largest_order+1) = 0;
    end
end

%% Return a grad structure
grad.info    = info;
grad.Gx      = Gx;    % [mT/m]
grad.Gy      = Gy;    % [mT/m]
grad.Gz      = Gz;    % [mT/m]
grad.R0      = R0;    % [m]
grad.lnorm   = lnorm; % [mT/m]
grad.alpha_z = alpha_z;
grad.alpha_x = alpha_x;
grad.alpha_y = alpha_y;
grad.beta_z  = beta_z;
grad.beta_x  = beta_x;
grad.beta_y  = beta_y;

end