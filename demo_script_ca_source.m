clear;
close all;
%% load file

%addpath(genpath('utilities'));
File=44;
filelist=dir('*.tif');
for File=1:length(filelist)

nam = filelist(File).name;
%nam = 'C:\Temp\bb\f2c1.tif';          % insert path to tiff stack here
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=3000;					% user input: how many frames to read   (optional, default until the end)

Y = bigread2(nam,sframe,num2read);
%Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to single

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels



%% Set parameters

K = 900;                                           % number of components to be found
tau =[3,3];                                          % std of gaussian kernel (size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.9;                                  % merging threshold
options = [];
options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'init_method','sparse_NMF',...                      % initialize algorithm with plain NMF  
    'snmf_max_iter',500,...                     % maximum number of iterations
    'search_method','dilate',...       % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold %'init_method','HALS',...    
    'tsub',1,...
    'ssub',1,...        %'conn_comp',true,...
    'conn_comp',false,...               %'spatial_method','constrained',...    
    'nb',2,...
    'noise_range',[0.25 0.75],...
    'spatial_method','constrained',...
    'rem_prct',5,...
    'df_window',40,...
    'medw',[2,2],...
    'gSig',tau...,
    );

%% Data pre-processing

[P,Y] = preprocess_data(Y,p);

%% fast initialization of spatial components using greedyROI and HALS
[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize

% display centers of found components
Cn =  correlation_image(Y); %reshape(P.sn,d1,d2);  %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
% figure;imagesc(Cn);
%     axis equal; axis tight; hold all;
%     scatter(center(:,2),center(:,1),'mo');
%     title('Center of ROIs found from initialization algorithm');
%     drawnow;

%% manually refine components (optional)
refine_components = false;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end
    
%% update spatial components
Yr = reshape(Y,d,T);
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

%[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

%% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

% %% classify components
% [ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(Y,A,C,b,f,YrA,options);
% 
% run_GUI = true;
% if run_GUI
%     Coor = plot_contours(A,Cn,options,1); close;
%     GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);   
%     options = GUIout{2};
%     keep = GUIout{3};    
% end

%% merge found components
[Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

% %%
% display_merging = 0; % flag for displaying merging example
% if and(display_merging, ~isempty(merged_ROIs))
%     i = 1; %randi(length(merged_ROIs));
%     ln = length(merged_ROIs{i});
%     figure;
%         set(gcf,'Position',[300,300,(ln+2)*300,300]);
%         for j = 1:ln
%             subplot(1,ln+2,j); imagesc(reshape(A(:,merged_ROIs{i}(j)),d1,d2)); 
%                 title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
%         end
%         subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
%                 title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
%         subplot(1,ln+2,ln+2);
%             plot(1:T,(diag(max(C(merged_ROIs{i},:),[],2))\C(merged_ROIs{i},:))'); 
%             hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
%             title('Temporal Components','fontsize',16,'fontweight','bold')
%         drawnow;
% end

%% repeat
Pm.p = p;    % restore AR value
[A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);

%% do some plotting

[A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
K_m = size(C_or,1);
[C_df,~] = extract_DF_F(Yr,A_or,C_or,P_or,options); % extract DF/F values (optional)

% %contour_threshold = 0.95;                       % amount of energy used for each component to construct contour plot
% figure;
% [Coor,json_file] = plot_contours(A_or,Cn,options,1); % contour plot of spatial footprints
% %savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)

%% display components

%plot_components_GUI(Yr,A_or,C_or,b2,f2,Cn,options)

%% make movie


Good_Components(File).DF=C_df;
Good_Components(File).RawCalcium=C_or;
Good_Components(File).name=filelist(File).name;
Good_Components(File).ROI=A_or;
Good_Components(File).ROIdims=[d1 d2];
Good_Components(File).P=P_or;
Good_Components(File).Spikes=S_or;
Good_Components(File).Background=f2;
Good_Components(File).Cn=Cn;
end



figure;
for i=1:length(filelist)
    subplot(6,4,i);imagesc(Good_Components(i).DF);colormap hot
end