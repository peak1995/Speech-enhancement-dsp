function [Output, NumSegments] = KFrame(Input, WindowLength, Window, HoppingSize)
% Chopper windows the signal based on window length, shift percantage and
% uses Hamming windowing technique.

% Number of samples to hop.
HoppingSamples = fix(WindowLength.*HoppingSize);

% Number of segments.
NumSegments = fix(((length(Input)-WindowLength)/HoppingSamples) + 1);

% Index matrix which guides the signal through chopping process.
Index = (repmat(1:WindowLength,NumSegments,1) + repmat((0:(NumSegments-1))'*HoppingSamples,1,WindowLength))';

% Final window which multiplies with original signal to give pieces of it.
FinalWindow = repmat(Window,1,NumSegments);

% Ta-da... 
Output = Input(Index).*FinalWindow;
end