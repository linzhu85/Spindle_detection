function D_fun_spindle_detection_wav(subject_id,run,Nap)
   % direct to the source file
   
   %Nap=0; % 0 indicates waking data while 1 indicates Nap data

   root_dir='/autofs/cluster/manoach/K24_MEG_EEG_Study/Analysis_Spindles';
   %subject_id='KC_01';run='02';
   subject_group=subject_id(1:2);
   
if Nap==0
   source_dir=strcat(root_dir,'/Waking_Spindles/Waking_matfiles','/',subject_group);
   cd(source_dir);
   load(strcat(subject_id,'_',run,'_ArtRej.mat'))
else
    source_dir=strcat(root_dir,'/Preliminary_Results_wavelet'); 
    cd(source_dir);
    load(strcat(subject_id,'_N2_MST_fast.mat'));
end

   % read the sampling rate 
   Fs=hdr.info.sfreq;
    
   % Time that the signal has to exceed background noise
    tthresh = .4 ; % sec, default is 400msec

   % Factor to multiply the mean of the signal in the spindle detection
    ampl_factor = 4.5; % default is 4.5    
    
   %%%%%%%wavelet setting%%%%%%%%%%%

   % Mother wavelet center frequency
    fc = 13.5;
    % Mother wavelet bandwidth
    fb = .5; 
    
    % Seconds before and after detection to determine spindle parameters
    segmDurBef = 2; % in secs
    segmDurAft = 2;

    % The appropriate scale depends on sampling frequency -- Use function 'SCAL2FRQ' to figure out the relationship between frequency and scale at your sampling rate
    % scale = 7.4; % Specifies scale to examine at 200Hz -- 
    
    % Feeding in the frequency gives you the scale, feeding in the scale gives you the frequency
    scale = scal2frq(13.5,['cmor' num2str(fc) '-' num2str(fb)],1/Fs);
%   freq = scal2frq(scale,['cmor' num2str(fc) '-' num2str(fb)],1/Fs);
    
    nchannel=size(data,1); % define number of channels
    
    EEG_index=find([hdr.info.chs.kind]==2); % only consider the EEGs so far
    Data=data(EEG_index,:); % get all EEG channels

    %% Calculate the continuous wavelet transform for each channel
    for x=1:length(EEG_index)
        current_data=Data(x,:);
        EEGWave=cwt(current_data,scale,['cmor' num2str(fc) '-' num2str(fb)]); % conducts the wavelet tranformation on workspace variable "EEG"
        EEGChannel = real(EEGWave.^2); % Takes only the real component of the coefficients
        EEGChannel=EEGChannel';
        Data_new(:,x) = EEGChannel;
    end

    %% Set the Sampling Frequency and take Moving Average
    Data_new=Data_new.^2;
    window = ones((Fs/10),1)/(Fs/10); % create 100ms window to convolve with
    Data2 = filter(window,1,Data_new); % take the moving average using the above window

    %% Compute Threshhold (Separately for Each Channel)
    signalmean=mean(Data2); % compute mean amplitude of rectified signal
    threshold = signalmean.*ampl_factor; % % defines the threshold
   
    for ch=1:length(EEG_index) % Loop for each signal
        
        fprintf('Working on Channel %s.\n',hdr.info.ch_names{EEG_index(ch)});
        
        %% Detect spindles
        current_data = Data2(:,ch);

        if any(isnan(current_data))
            spindle_det(ch).bads = 1;
        else
            spindle_det(ch).bads = 0;
        end
        
        over = current_data>threshold(1,ch); % Mark all points over threshold as '1'
        locs = (zeros(1,length(current_data)))'; % Create a vector of zeros the length of the MS signal

        for i=1:((length(current_data))-(Fs*tthresh)); % for the length of the signal, if the sum of 40 concurrent points = Fs*0.4, mark a spindle
            if sum(over(i:(i+((Fs*tthresh)-1))))==(Fs*tthresh);
                locs(i,1)=1;
            end
        end

        spin=zeros((length(locs)),1);  % only mark a spindle in vector 'spin' at the end of a 400ms duration peak
        for i=1:length(locs);
            if locs(i,1)==1 && locs(i+1,1)==0;
                spin(i,1)=1;
            end
        end
        
        % added 9/17/2012: for every spindle marked in 'spin', delete the spindle if there is also a spindle within the second preceeding it.
        for i=(Fs+1):length(spin);  
            if spin(i,1)==1 && sum(spin((i-Fs):(i-1)))>0;
                spin(i,1)=0;
            end
        end

        spindle_det(ch).label = hdr.info.ch_names{ch};
        spindle_det(ch).sample = find(spin);
        spindle_det(ch).spindle_count = sum(spin);
        spindle_det(ch).backgr_mean = signalmean(ch);
        
        %% Discard spindles that are closer than 2 sec to the boundaries
        bound_idx = find(spindle_det(ch).sample + Fs*segmDurAft > size(Data,2) | spindle_det(ch).sample - Fs*segmDurAft < 1);
        if ~isempty(bound_idx)
            spindle_det(ch).sample(bound_idx) = [];
            spindle_det(ch).spindle_count = length(spindle_det(ch).sample); 
        end
        
        %% Calculate Spindle Parameters
        numSpindles = spindle_det(ch).spindle_count;
        numSpindles(isnan(numSpindles)) = 0;
%         
%         if numSpindles == 0 % no spindles on this channel
%             continue;
%         else

            segms = cell(numSpindles,1);
            % Make 4 sec segments around each spindle event
            for jj = 1:numSpindles
                startSmp = max(spindle_det(ch).sample(jj)-Fs*segmDurBef+1,1);
                endSmp = min(spindle_det(ch).sample(jj)+Fs*segmDurAft,length(Data));
                segms{jj} = Data(ch,startSmp:endSmp);
            end
            
            % get rid of truncated segments
            tmp = cellfun(@length, segms); % Find shorter segments
            segms(tmp ~= Fs*(segmDurBef+segmDurAft)) = {NaN*ones(1,Fs*(segmDurBef+segmDurAft))}; %Turn their data to nans 
            segms = cell2mat(segms);
            
            %% Load the results for this record
             if Nap==0
              cd /autofs/cluster/manoach/K24_MEG_EEG_Study/Analysis_Spindles/Waking_Spindles/Wavelet_Results
              file_name=strcat(subject_id,'_',run,'_MST_fast.mat'); % find the corresponding output
              load(file_name);
          else
              cd /autofs/cluster/manoach/K24_MEG_EEG_Study/Analysis_Spindles/Preliminary_Results_wavelet
              file_name=strcat(subject_id,'_N2_MST_fast.mat');
          end
            spindle_det=allSpindles;
             %% Discard spindles that overlap
            if length(spindle_det(EEG_index(ch)).startSample)>=2
                overl_idx = find(spindle_det(ch).startSample(2:end)-spindle_det(ch).endSample(1:end-1)<0);
                
                spindle_det(ch).sample(overl_idx)      = [];
                spindle_det(ch).energy(overl_idx)      = [];
                spindle_det(ch).peakAmp(overl_idx)     = [];
                spindle_det(ch).peakLoc(overl_idx)     = [];
                spindle_det(ch).startSample(overl_idx) = [];
                spindle_det(ch).endSample(overl_idx)   = [];
                spindle_det(ch).duration(overl_idx)    = [];
                spindle_det(ch).spindle_count = length(spindle_det(ch).duration);
          % find the original result and add one column for this parameter
            %%added 04/04/2018, get peak frequency of each spindles segms
            foo=[]; peakloc={};peakfreq={};
            for k=1:size(segms,1)
                foo=segms(k,:);
                freq_fft=abs(fft(foo));% get the FFT of spindle seg
                % needs to find the peak points between 12 to 15Hz  
                [m,n]=max(freq_fft(1,48:60));% freq resolution is Fs/N =0.25, 
                peakloc{k}=(48+n-2)*0.25; % all the spindle peak frquency saved as cell
                %f=Fs*(0:length(foo)-1)/length(foo);
                %plot(f(1:length(foo)/4),freq_fft(1:length(foo)/4));
            end
            peakfreq{ch}=peakloc; 
                     
%            
          
          
          if Nap==0
              cd /autofs/cluster/manoach/K24_MEG_EEG_Study/Analysis_Spindles/Waking_Spindles/Wavelet_Results
              file_name=strcat(subject_id,'_',run,'_MST_fast.mat'); % find the corresponding output
              load(file_name);
              allSpindles(EEG_index(ch)).peakfreq=peakfreq{ch};
              dest_name=strcat(subject_id,'_',run,'_MST_fast.mat'); % find the corresponding output
              save(dest_name,'allSpindles','param')
          else
              cd /autofs/cluster/manoach/K24_MEG_EEG_Study/Analysis_Spindles/Preliminary_Results_wavelet
              file_name=strcat(subject_id,'_N2_MST_fast.mat');
              allSpindles(EEG_index(ch)).peakfreq=peakfreq{ch};
              dest_name=strcat(subject_id,'_N2_MST_fast.mat');
              save(dest_name,'allSpindles','param')
          end
    end

%            %%%%%%%%%% get allspindles parameters%%%%%%%
%             [spindle_det(ch).duration, spindle_det(ch).energy, spindle_det(ch).peakAmp, spindle_det(ch).peakLoc, spindle_det(ch).startSample, spindle_det(ch).endSample] = spindle_duration_Energy_K24_vOrig_Fs200(segms, Fs); % ******
%             
%             % Convert to whole-signal timeline (samples given relative to 1 sec before detection point)
%             spindle_det(ch).startSample  = spindle_det(ch).startSample + spindle_det(ch).sample' - Fs*(segmDurBef-1)*ones(size([spindle_det(ch).startSample]));
%             spindle_det(ch).endSample    = spindle_det(ch).endSample   + spindle_det(ch).sample' - Fs*(segmDurBef-1)*ones(size([spindle_det(ch).endSample]));
%             spindle_det(ch).peakLoc      = spindle_det(ch).peakLoc     + spindle_det(ch).sample' - Fs*(segmDurBef-1)*ones(size([spindle_det(ch).peakLoc]));
%             
%             %% Discard spindles that long less than x samples
%             spindle_det(ch).sample([spindle_det(ch).duration]<tthresh)      = [];
%             spindle_det(ch).energy([spindle_det(ch).duration]<tthresh)      = [];
%             spindle_det(ch).peakAmp([spindle_det(ch).duration]<tthresh)     = [];
%             spindle_det(ch).peakLoc([spindle_det(ch).duration]<tthresh)     = [];
%             spindle_det(ch).startSample([spindle_det(ch).duration]<tthresh) = [];
%             spindle_det(ch).endSample([spindle_det(ch).duration]<tthresh)   = [];
%             spindle_det(ch).duration([spindle_det(ch).duration]<tthresh)    = [];
%             spindle_det(ch).spindle_count = length(spindle_det(ch).duration);
%             
%             %% Discard spindles that overlap
%             if length(spindle_det(ch).startSample)>=2
%                 overl_idx = find(spindle_det(ch).startSample(2:end)-spindle_det(ch).endSample(1:end-1)<0);
%                 
%                 spindle_det(ch).sample(overl_idx)      = [];
%                 spindle_det(ch).energy(overl_idx)      = [];
%                 spindle_det(ch).peakAmp(overl_idx)     = [];
%                 spindle_det(ch).peakLoc(overl_idx)     = [];
%                 spindle_det(ch).startSample(overl_idx) = [];
%                 spindle_det(ch).endSample(overl_idx)   = [];
%                 spindle_det(ch).duration(overl_idx)    = [];
%                 spindle_det(ch).spindle_count = length(spindle_det(ch).duration);
                
          
      
    %end % End the loop for each channel
