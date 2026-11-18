function plot_streamlines_from_nucleus(Nuc2Plot,Color2Use,AlphaData)

% -------------------------------------------------------
% Plot all streamlines from a given thalamic nucleus
%
% INPUTS:
% tracks        - mrtrix tracks structure (from read_mrtrix_tracks)
% Assignments   - [N x 2] matrix (thal origin, cortical end)
% nucleusID     - nucleus index (after remapping, 1–35)
% Colors        - Nx3 colour matrix (0–1)
% Thal2Ignore   - nuclei to ignore (vector)
% -------------------------------------------------------

tracks = read_mrtrix_tracks ('./ExampleIndvTract/combined_edited_tracts_excl_subMask.tck');
track_assignments = dlmread('./ExampleIndvTract/CC01236XX16_155830_thalEdits_dwi_rand100_assignments.txt');

tckValid = find(track_assignments(:,1)==(Nuc2Plot+200) & track_assignments(:,2)>100 & track_assignments(:,2)<=200);

%axis tight

ROI = track_assignments(tckValid,2)-100;

Ntracks = length(tckValid);

tracks_data = cell(Ntracks,1);

totalTractLengths = 0;

for i = 1:Ntracks
    tracks_data{i} = tracks.data{tckValid(i)};
totalTractLengths = totalTractLengths + 1 + size(tracks_data{i},1);

end

% Initialize arrays to store streamlines and colors
STREAMLINES = zeros(totalTractLengths, 3);
STREAMLINES_COLOR = zeros(totalTractLengths, 3);
STREAMLINES_ALPHADATA = zeros(totalTractLengths, 1);
STREAMLINES_IND = zeros(totalTractLengths, 1);
STREAMLINES_DIR = zeros(totalTractLengths, 3);
% Create another waitbar to track progress
f = waitbar(0, 'Starting');
IND = 0;
TRACT_ID = 1;
%turbo_cmap = turbo(250);

StepRate = 1;
    
    % Loop through tracks associated with a seed
    for i = 1:Ntracks
        tract = tracks_data{i};
        TractLength = length(tract) + IND;
        %color2use = Colors(ROIind(ROI(i)),:);
        %color2use = turbo_cmap(ColorIND(j),:);
        
        % Append streamline coordinates and colors
        STREAMLINES(IND+1:TractLength+1, :) = [[tract(:, 1)' NaN]', [tract(:, 2)' NaN]', [tract(:, 3)' NaN]'];
        STREAMLINES_COLOR(IND+1:TractLength+1, :) = repmat(Color2Use, size(tract, 1) + 1, 1);

        %STREAMLINES_ROI(IND+1:TractLength+1, :) = ROI(i);
        STREAMLINES_IND(IND+1:TractLength+1, :) = TRACT_ID;
        STREAMLINES_ALPHADATA(IND+1:TractLength+1, 1) = AlphaData(ROI(i));

        % joint = ([tract(2:end,:); NaN NaN NaN]-tract);
        % joint(end,:) = joint(end-1,:);
        % temp_joint = joint;
        % joint(:,1) = temp_joint(:,1);
        % joint(:,2) = temp_joint(:,2);
        % cdata = [abs(joint./StepRate); NaN NaN NaN];
        % STREAMLINES_DIR(IND+1:TractLength+1,:) = cdata;


        % Update index for next tract
        IND = TractLength + 1;
        TRACT_ID = TRACT_ID + 1;
        % Update waitbar progress
    waitbar(i / Ntracks, f, ['Finished tract ', num2str(i), ' of ', num2str(Ntracks)]);
    end
    
    

close(f)

% Create a figure for tractogram visualization
%figure('Position', [0 0 1920 1080])

%figure

cdata_all = reshape(STREAMLINES_COLOR, length(STREAMLINES_COLOR), 1, 3);   

%cdata_all = reshape(STREAMLINES_DIR, length(STREAMLINES_DIR), 1, 3);   

streamlines = patch(STREAMLINES(:, 1), STREAMLINES(:, 2), STREAMLINES(:, 3), 0);
set(streamlines, 'CData', cdata_all, 'EdgeColor', 'interp', 'FaceColor', 'interp', 'Clipping', 'off') 

% Configure lighting and view settings for 3D visualization
camlight(80, -10);     % Add light from upper-right direction
camlight(-80, -10);    % Add light from upper-left direction
view([90 0]);          % Set the initial view to a lateral view

% Configure axis settings for visualization
axis off;              % Turn off axis display
axis tight;            % Set the axis limits to fit the data
axis equal;            % Set equal scaling for all axis
axis vis3d;            % Set 3D visualization mode

%streamlines.EdgeAlpha=.15;

streamlines.EdgeAlpha='flat';
streamlines.FaceVertexAlphaData=STREAMLINES_ALPHADATA;
streamlines.AlphaDataMapping='none';
hold on

RIGHT = gifti('./ExampleIndvTract/week-40_hemi-right_space-sub-CC01236XX16_ses-155830_dwi_wm.surf.gii');

surfaceR.vertices = RIGHT.vertices;
surfaceR.faces = RIGHT.faces;

pMNIsurfaceR = patch(surfaceR);
set(pMNIsurfaceR,'EdgeColor','none','FaceColor',[.5 .5 .5],'Clipping','off');
pMNIsurfaceR.FaceLighting = 'gouraud';
material dull
pMNIsurfaceR.FaceAlpha = 0.1;
