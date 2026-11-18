load('DemographicDataQCpass.mat')

% --- Age range in weeks/days ---
minAge = min(DemographicDataQCpass.ScanAge);
maxAge = max(DemographicDataQCpass.ScanAge);

% Convert fractional weeks to [weeks, days]
ScanAgerangeStr = sprintf('%.0f weeks %.0f days – %.0f weeks %.0f days', ...
    floor(minAge), round((minAge - floor(minAge))*7), ...
    floor(maxAge), round((maxAge - floor(maxAge))*7));

% --- Median scan age ---
medAge = median(DemographicDataQCpass.ScanAge);
ScanAgemedStr = sprintf('%.0f weeks %.0f days', ...
    floor(medAge), round((medAge - floor(medAge))*7));

% --- Age range (birth) in weeks/days ---
minAge = min(DemographicDataQCpass.BirthAge);
maxAge = max(DemographicDataQCpass.BirthAge);

% Convert fractional weeks to [weeks, days]
BirthAgerangeStr = sprintf('%.0f weeks %.0f days – %.0f weeks %.0f days', ...
    floor(minAge), round((minAge - floor(minAge))*7), ...
    floor(maxAge), round((maxAge - floor(maxAge))*7));

% --- Median birth age ---
medAge = median(DemographicDataQCpass.BirthAge);
BirthAgemedStr = sprintf('%.0f weeks %.0f days', ...
    floor(medAge), round((medAge - floor(medAge))*7));

% --- Count females ---
if iscellstr(DemographicDataQCpass.sex) || isstring(DemographicDataQCpass.sex)
    nFemales = sum(strcmpi(DemographicDataQCpass.sex,'female'));
elseif iscategorical(DemographicDataQCpass.sex)
    nFemales = sum(DemographicDataQCpass.sex == 'female');
else
    error('Sex variable should be string or categorical');
end

% --- Print results ---

fprintf('Birth age range: %s\n', BirthAgerangeStr);
fprintf('Median birth age: %s\n', BirthAgemedStr);

fprintf('Scan age range: %s\n', ScanAgerangeStr);
fprintf('Median scan age: %s\n', ScanAgemedStr);
fprintf('Total number of females: %d\n', nFemales);