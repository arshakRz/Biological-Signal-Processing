%% Q3 Part one: segmentation using classification
%%

PCGCellArray = train_recordings;
annotationsArray = train_annotations;
numberOfStates = 4;
numPCGs = length(PCGCellArray);
Fs = 1000;
% A matrix of the values from each state in each of the PCG recordings:
state_observation_values = cell(numPCGs,numberOfStates);


for PCGi = 1:length(PCGCellArray)
    PCG_audio = PCGCellArray{PCGi};
    
    S1_locations = annotationsArray{PCGi,1};
    S2_locations = annotationsArray{PCGi,2};
    
    [PCG_Features, featuresFs] = getSpringerPCGFeatures(PCG_audio, Fs);
    
    PCG_states = labelPCGStates(PCG_Features(:,1),S1_locations, S2_locations, featuresFs, 1);
    %% your code here:



    %%
end
    
    