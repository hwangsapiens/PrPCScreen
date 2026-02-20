% MATLAB script to plot a Manhattan plot from data in an Excel file

% Define the Excel file and sheet name
filename = '\\fs-home\wanha$\Desktop\PrPcScreenFullResults\PrP_genes_and_NT_ordered_with_chromosome.xlsx';  % Replace with the path to your file
sheetname = 'manhattanplot2';

% Read data from the Excel file
data = readtable(filename, 'Sheet', sheetname);

% Extract columns
geneSymbols = data.Gene_symbol;
log2FC = data.Mean_log2FC;
chromosomes = data.Chromosome; % Chromosome labels in the format 'chr1', 'chr2', ...
positions = data.Start_Position;

% Convert chromosome names to numeric for plotting
chromosomes_numeric = zeros(size(chromosomes));
uniqueChromosomes = unique(chromosomes, 'stable');

for i = 1:length(uniqueChromosomes)
    if strcmp(uniqueChromosomes{i}, 'ChrX')
        chromosomes_numeric(strcmp(chromosomes, uniqueChromosomes{i})) = 23;  % X chromosome
    elseif strcmp(uniqueChromosomes{i}, 'ChrY')
        chromosomes_numeric(strcmp(chromosomes, uniqueChromosomes{i})) = 24;  % Y chromosome
    else
        % Extract numeric part for autosomes (chr1, chr2, ..., chr22)
        numStr = regexprep(uniqueChromosomes{i}, 'Chr', ''); % Remove 'Chr'
        chromosomes_numeric(strcmp(chromosomes, uniqueChromosomes{i})) = str2double(numStr);
    end
end

% Calculate cumulative position for each chromosome to spread chromosomes along x-axis
uniqueChromosomeNumbers = unique(chromosomes_numeric);
chromOffsets = zeros(length(uniqueChromosomeNumbers), 1);

% Calculate offsets and adjusted positions
adjusted_positions = zeros(size(positions));

for i = 1:length(uniqueChromosomeNumbers)
    chr = uniqueChromosomeNumbers(i);
    indices = find(chromosomes_numeric == chr);
    
    if ~isempty(indices)
        % Get positions for this chromosome
        chr_positions = positions(indices);
        
        % If it's not the first chromosome, update offset
        if i > 1
            chromOffsets(i) = chromOffsets(i-1) + max(positions(chromosomes_numeric == uniqueChromosomeNumbers(i-1)));
        end
        
        % Adjust the positions for the current chromosome
        adjusted_positions(indices) = chr_positions + chromOffsets(i);
    end
end

% Plot the Manhattan plot
figure;
hold on;

% Define colors using hex codes
colors = [hex2dec('9a')/255, hex2dec('7e')/255, hex2dec('6f')/255; 
          hex2dec('54')/255, hex2dec('47')/255, hex2dec('3f')/255];  % Colors from hex codes

% Plot each chromosome with the defined colors
for i = 1:length(uniqueChromosomeNumbers)
    chrIndices = chromosomes_numeric == uniqueChromosomeNumbers(i);
    scatter(adjusted_positions(chrIndices), log2FC(chrIndices), 10, ...
        'MarkerEdgeColor', colors(mod(i-1,2)+1, :), 'MarkerFaceColor', colors(mod(i-1,2)+1, :));
end

% Highlight the PRNP gene
prnpIndex = find(strcmp(geneSymbols, 'PRNP')); % Find the index for PRNP gene
if ~isempty(prnpIndex)
    prnpColor = [hex2dec('BC')/255, hex2dec('7C')/255, hex2dec('7C')/255]; % Convert hex to RGB
    scatter(adjusted_positions(prnpIndex), log2FC(prnpIndex), 50, ...
        'MarkerEdgeColor', prnpColor, 'MarkerFaceColor', prnpColor, 'LineWidth', 1.5); % Highlight PRNP
    text(adjusted_positions(prnpIndex), log2FC(prnpIndex), ' PRNP', ...
        'FontSize', 10, 'FontName', 'Arial', 'Color', 'black', 'VerticalAlignment', 'bottom'); % Label PRNP
end

% Add labels and format the plot
xlabel('Chromosome', 'FontName', 'Arial', 'FontSize', 12); % Set Arial font
ylabel('Mean log2 Fold Change', 'FontName', 'Arial', 'FontSize', 12); % Set Arial font
title('Manhattan Plot', 'FontName', 'Arial', 'FontSize', 12); % Set Arial font

% Set tick positions and labels
xticks(chromOffsets(1:end-1) + diff(chromOffsets) / 2);
% Create labels as 1, 2, ..., X, Y
xticklabels(arrayfun(@(x) x, 1:22, 'UniformOutput', false)); % Chromosomes 1-22
xticklabels{23} = 'X'; % Set label for chromosome 23
xticklabels{24} = 'Y'; % Set label for chromosome 24
ax = gca; % Get current axes
ax.XTickLabel = str2double(ax.XTickLabel); % Convert XTickLabels to numeric format
set(gca, 'FontName', 'Arial', 'FontSize', 10); % Set font and size for the ticks

% Adjust the plot limits for better visualization
ylim([min(log2FC) - 1, max(log2FC) + 1]);
xlim([0, max(adjusted_positions)]);

hold off;
