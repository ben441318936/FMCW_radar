%% load/import raw data
data = 'rawData3D_simple2D'; % Change only this line
raw_data = load(data);
raw_data = raw_data.(data);

%% declare some vals we will use.
nFFT=1024; %space and time fft n are both 1024 so just package into 1 var
targetZ=280; %range of target, measure in z axis as a range for img slice 

%67.5 cm depth - 1 cm

sampleX=200/406;
sampleY=2;
light_spd=physconst('lightspeed');
fs=9121e3; %sample rate
ts=1/fs; %sample period
k=63.343e12; %slope const in (hz/sec)  ??

%% range-fft rawdata

raw_data_fft=fft(raw_data,nFFT);

%% Range focusing to z0
I_delay = 4.5225e-10; % Instrument delay for range calibration (corresponds to a 6.78cm range offset)
range_bin = round(k*ts*(2*targetZ*10^-3/light_spd+I_delay)*nFFT); % corresponing range bin
sar_data = squeeze(raw_data_fft(range_bin+1,:,:));

%% do match filter 
f0=77e9; %initial frequency
x=sampleX*(-(nFFT-1)/2:(nFFT-1)/2)*1e-3;%xsample*[-512:512] into mm
y=(sampleY*(-(nFFT-1)/2:(nFFT-1)/2)*1e-3).';%ysample*[-512:512] into mm but transposed bc vertical axis
z=targetZ*1e-3;
k1=2*pi*f0/light_spd;
matched_filter=exp(-1i*2*k1*sqrt(x.^2+y.^2+z^2));


%% reconstruct sar image
im_size = 200; % Size of image area in mm

%y x x
[y_sar_f,x_sar_f]=size(sar_data);
[y_match_f,x_match_f]=size(matched_filter);

%equalize dims w/ zero padding so consistent w/ each other
if(x_match_f>x_sar_f)
    sar_data=padarray(sar_data,[0 floor((x_match_f-x_sar_f)/2)],0,'pre');
    sar_data=padarray(sar_data,[0 ceil((x_match_f-x_sar_f)/2)],0,'post');
else
    matched_filter=padarray(matched_filter,[0 floor((x_sar_f-x_match_f)/2)],0,'pre');
    matched_filter=padarray(matched_filter,[0 ceil((x_sar_f-x_match_f)/2)],0,'post');
end

if(y_match_f>y_sar_f)
    sar_data=padarray(sar_data,[floor((y_match_f-y_sar_f)/2) 0],0,'pre');
    sar_data=padarray(sar_data,[ceil((y_match_f-y_sar_f)/2) 0],0,'post');
else
    matched_filter=padarray(matched_filter,[floor((y_sar_f-y_match_f)/2) 0],0,'pre');
    matched_filter=padarray(matched_filter,[ceil((y_sar_f-y_match_f)/2) 0],0,'post');
end

%apply FFTs and create the actual SAR IMAGE
sar_data_FFT=fft2(sar_data);
matched_filter_FFT=fft2(matched_filter);
final_sar_image=fftshift(ifft2(sar_data_FFT .* matched_filter_FFT));

[y_target_t,x_target_t]=size(final_sar_image);

X=sampleX* (-(x_target_t-1)/2: (x_target_t-1)/2); %base on size of im but is also in mm alrdy
Y=sampleX* (-(x_target_t-1)/2: (x_target_t-1)/2);

%crop to img of related region
indX=X>(-im_size/2)&X<(im_size/2); 
indY=Y>(-im_size/2)&Y<(im_size/2);

X=X(indX);
Y=Y(indY);
final_sar_image=final_sar_image(indY,indX);

%finally, plot sar image.
figure;
mesh(X,Y,abs(fliplr(final_sar_image)),'FaceColor','interp','LineStyle','none')
view(2)
colormap('jet')
xlabel('Horizontal (mm)')
ylabel('Vertical (mm)')
figure_title = "SAR Image - Matched Filter Response";
title(figure_title)


