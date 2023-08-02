%% 1. Simulate sparse basis functions of the natural images:
clc; clear; close all;
addpath 'sparsenet/'
load IMAGES.mat

num_trials = 3000;
A = rand(256, 144) - 0.5;
A = A * diag(1./sqrt(sum(A.*A)));
figure(1)
colormap(gray)
sparsenet

figure(1)
set(gcf, 'PaperPositionMode', 'auto')
print("BasisFunc_" + size(A, 1) + "_" + size(A, 2), '-dpng', '-r0')

figure(2)
set(gcf, 'PaperPositionMode', 'auto')
print("NormBasisFunc_" + size(A, 1) + "_" + size(A, 2), '-dpng', '-r0')

%% 2. Study the effect of different datasets:
%% Yale
clc; clear; close all;
addpath 'sparsenet/'

IMAGES = zeros(192, 168, 10);

for i = 1:10
    if i < 10
        img = imread("CroppedYalePNG/yaleB0" + num2str(i) + "_P00A+000E+20.png");
    else
        img = imread("CroppedYalePNG/yaleB" + num2str(i) + "_P00A+000E+20.png");
    end
    IMAGES(:,:,i) = img;
end

% Whitening images
for img_index = 1:size(IMAGES, 3)
    img = IMAGES(:,:,img_index);
    img = img/max(img,[],'all');  
    img = img-mean(img, 'all');
    IMAGES_tmp(:,:,img_index) = img(13:end-12,:);
end
IMAGES = IMAGES_tmp;
IMAGES = whiteningFunc(IMAGES);

num_trials = 3000;
A = rand(64, 144) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));
figure(1)
colormap(gray)
sparsenet

figure(1)
set(gcf, 'PaperPositionMode', 'auto')
print("YaleBasisFunc_" + size(A, 1) + "_" + size(A, 2), '-dpng', '-r0')

figure(2)
set(gcf, 'PaperPositionMode', 'auto')
print("YaleNormBasisFunc_" + size(A, 1) + "_" + size(A, 2), '-dpng', '-r0')

%% MNIST
clc; clear; close all;
addpath 'sparsenet/'

IMAGES = load('mnist-original.mat');
IMAGES.data = reshape(IMAGES.data, [28, 28, 70000]);
IMAGES = IMAGES.data(:,:,1:10);

% Whitening images
for img_index = 1:size(IMAGES, 3)
    img = IMAGES(:,:,img_index);
    img = img/max(img,[],'all');  
    img = img-mean(img, 'all');
    IMAGES(:,:,img_index) = img;
end
IMAGES = whiteningFunc(IMAGES);

num_trials = 3000;
A = rand(64, 144) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));
figure(1)
colormap(gray)
sparsenet

figure(1)
set(gcf, 'PaperPositionMode', 'auto')
print("MnistBasisFunc_" + size(A, 1) + "_" + size(A, 2), '-dpng', '-r0')

figure(2)
set(gcf, 'PaperPositionMode', 'auto')
print("MnistNormBasisFunc_" + size(A, 1) + "_" + size(A, 2), '-dpng', '-r0')

%% caltech 101
clc; clear; close all;
addpath 'sparsenet/'

IMAGES = zeros(200, 200, 10);

for i = 1:10
    if i < 10
        img = imread("caltech-101/ant/image_000" + num2str(i) + ".jpg");
    else
        img = imread("caltech-101/ant/image_00" + num2str(i) + ".jpg");
    end
    if size(img, 3) > 1
        img = rgb2gray(img);
    end
    img = imresize(img, [200, 200]);
    IMAGES(:,:,i) = img;
end

% Whitening images
for img_index = 1:size(IMAGES, 3)
    img = IMAGES(:,:,img_index);
    img = img/max(img,[],'all');  
    img = img-mean(img, 'all');
    IMAGES(:,:,img_index) = img;
end
IMAGES = whiteningFunc(IMAGES);

num_trials = 3000;
A = rand(64, 144) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));
figure(1)
colormap(gray)
sparsenet

figure(1)
set(gcf, 'PaperPositionMode', 'auto')
print("CaltechBasisFunc_" + size(A, 1) + "_" + size(A, 2), '-dpng', '-r0')

figure(2)
set(gcf, 'PaperPositionMode', 'auto')
print("CaltechNormBasisFunc_" + size(A, 1) + "_" + size(A, 2), '-dpng', '-r0')

%% 3. Study the dynamics of the sparse coefficients
clc; clear; close all;
addpath 'sparsenet/'

v = VideoReader('BIRD.avi');
v_frames_num = v.Duration * v.FrameRate;

channels_num = 1; % Grayscale frames
frame = zeros(v.Height, v.Width, channels_num, v_frames_num);
frame_rgb = zeros(v.Height, v.Width,3, v_frames_num);

im_size = size(frame, 1);

k = 1;
while hasFrame(v)
    frame_rgb(:,:,:,k) = readFrame(v);
    frame(:,:,:,k) = im2gray(uint8(frame_rgb(:,:,:,k))); 
    k = k+1;
end

frame = squeeze(frame);

% resize images
frame = frame(1:im_size, 1:im_size, :);
rgb_vid = frame_rgb(1:im_size, 1:im_size, :, :);
rgb_vid = uint8(rgb_vid);

frame = whiteningFunc(frame);
video_patches = frame;
 
IMAGES = video_patches(:,:,11:20);

num_trials = 1000;
batch_size = 100;

num_images=size(IMAGES,3);
image_size=size(IMAGES,1);
BUFF = 4;

A = rand(576,100) - 0.5;
A = A*diag(1./sqrt(sum(A.*A)));

[L M] = size(A);
sz = sqrt(L);

eta = 1.0;
noise_var = 0.01;
beta = 2.2;
sigma =0.316;
tol =.01;

VAR_GOAL = 0.1;
S_var = VAR_GOAL*ones(M,1);
var_eta = .001;
alpha = .02;
gain = sqrt(sum(A.*A));

X = zeros(L,batch_size);

display_every = 10;

h = display_network(A,S_var);

for t = 1:num_trials
    i = ceil(num_images*rand);
    this_image = IMAGES(:,:,i);
    if(length(find(isnan(this_image) == 1)) == 0)
        for i = 1:batch_size
            r = BUFF+ceil((image_size-sz-2*BUFF)*rand);
            c = BUFF+ceil((image_size-sz-2*BUFF)*rand);
            X(:,i) = reshape(this_image(r:r+sz-1,c:c+sz-1),L,1);
        end

        S = cgf_fitS(A,X,noise_var,beta,sigma,tol);
        E = X-A*S;

        dA = zeros(L,M);

        for i = 1:batch_size
            dA = dA + E(:,i)*S(:,i)';
        end

        dA = dA/batch_size;
        A = A + eta*dA;

        for i=1:batch_size
            S_var = (1-var_eta)*S_var + var_eta*S(:,i).*S(:,i);
        end
        gain = gain .* ((S_var/VAR_GOAL).^alpha);
        normA = sqrt(sum(A.*A));
        for i = 1:M
            A(:,i) = gain(i)*A(:,i)/normA(i);
        end

        if (mod(t,display_every) == 0)
            display_network(A,S_var,h);
        end
    end
end

figure(1)
set(gcf, 'PaperPositionMode', 'auto')
print("BirdBasisFunc_" + size(A, 1) + "_" + size(A, 2), '-dpng', '-r0')

figure(2)
set(gcf, 'PaperPositionMode', 'auto')
print("BirdNormBasisFunc_" + size(A, 1) + "_" + size(A, 2), '-dpng', '-r0')

%% S Bird
S_BIRD = zeros(size(video_patches,3)-10+1, size(S,1), size(S,2));
S_BIRD(1,:,:) = S;

beta = 2.2;
sigma = 0.316;
tol = .01;
for I = 11:size(video_patches,3)
    IMAGES = video_patches(:,:,I);

    num_images = size(IMAGES,3);
    image_size = size(IMAGES,1);

    this_image = IMAGES;
    X = zeros(size(A,1),(size(this_image,1)/sqrt(size(A,1)))^2);
    for i = 1:size(X,2)
        X(:,i) = reshape(this_image((int32((i*sqrt(size(A,1)))/image_size)+1):...
            (int32((i*sqrt(size(A,1)))/image_size)+24),...
            mod(i*sqrt(size(A,1)),image_size)+1:...
            mod(i*sqrt(size(A,1)),image_size)+24),...
            [sqrt(size(A,1))^2 1]);
    end

    S = cgf_fitS(A,X,noise_var,beta,sigma,tol);
    
    S_BIRD(I-9,:,:) = S;
end

%% Write videos
frames = [];

figure();
for i = 1:size(S_BIRD,1)
    axis square
    pcolor(1:144,1:100,squeeze(S_BIRD(i,:,:)));
    xlabel('WindowNum','interpreter','latex');
    ylabel('BasisFunc num','interpreter','latex');
    colorbar
    colormap gray
    caxis([min(S_BIRD(:)) max(S_BIRD(:))])
    frames = [frames getframe(gcf)];
    drawnow
end

writerObj = VideoWriter("Heatmap");
writerObj.FrameRate = 20;
writerObj.Quality = 100;

open(writerObj);
for i=1:length(frames)
    frame = frames(i) ;
    writeVideo(writerObj,frame);
end
close(writerObj)

frames = [];
figure();
selected_window = 60;
for i = 1:size(S_BIRD,1)
    plot(1:100,S_BIRD(i,:,selected_window),'LineWidth',1,'Color','b')
    ylim([-1 1])
    xlim([0 100])
    hold off
    xlabel('BasisFunc num','interpreter','latex');
    ylabel('Basis coeff val','interpreter','latex');
    grid on;
    frames = [frames getframe(gcf)];
    drawnow
end

writerObj = VideoWriter("CoeffVsBasisNum");
writerObj.FrameRate = 20;
writerObj.Quality = 100;

open(writerObj);
for i=1:length(frames)
    frame = frames(i) ;
    writeVideo(writerObj,frame);
end
close(writerObj)

function IMAGES = whiteningFunc(IMAGES)
    N = size(IMAGES, 1);
    M = size(IMAGES, 3);

    [fx fy]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);
    rho=sqrt(fx.*fx+fy.*fy);
    f_0=0.4*N;
    filt=rho.*exp(-(rho/f_0).^4);

    for i=1:M
        image = IMAGES(:,:,i);
        If=fft2(image);
        imagew = real(ifft2(If.*fftshift(filt)));
        IMAGES_tmp(:,i)=reshape(imagew,N^2,1);
    end

    IMAGES_tmp = sqrt(0.1) * IMAGES_tmp/sqrt(mean(var(IMAGES_tmp)));
    IMAGES = reshape(IMAGES_tmp, size(IMAGES));
end
