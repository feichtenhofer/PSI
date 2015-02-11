function s = PSI( I, percentile )
% PSI(I) determines the perceptual sharpness s of an image I
%
%Reference:
%   
%   C Feichtenhofer, H Fassold, P Schallauer
%   "A perceptual image sharpness metric based on local edge gradient
%   analysis", IEEE Signal Processing Letters, 20 (4), 379-382, 2013
%   
%
%   Written by Christoph Feichtenhofer (cfeichtenhofer AT gmail.com)
%   feichtenhofer.github.io   
%


%% Parameters
if nargin < 2
    percentile=22; % percentage of blocks to use for metric 
end
% BLOCKSIZE [Def:32] Size for averaging of width measurements
% THRESHOLD_W [Def:2] Sum of widths in block to proccess block further
BLOCKSIZE = 32;
THRESHOLD_W = 2;


if ( size(I,3) > 1 ) 
    I = rgb2gray(I);    
end

% sobel_tr [Def: []] If value is assigned, this is the constant sobel threshold,
sobel_tr = [];
[edges] = edge(I,'Sobel',sobel_tr);

I = double(I) / 255;
      
QUOTA_W = percentile/100; 

row_blocks = floor(size(I,1)/BLOCKSIZE);
col_blocks = floor(size(I,2)/BLOCKSIZE);


%% calculate angles and round them, then calc. horz/vert widths.
[m, n] = size(I);
edge_widths = zeros(m,n);
widths_count = 0;

% calculate gradient 
Ix = [I(:,2)-I(:,1), 0.5*(I(:,3:end)-I(:,1:end-2)), I(:,end)-I(:,end-1)];
Iy = [I(2,:)-I(1,:); 0.5*(I(3:end,:)-I(1:end-2,:)); I(end,:)-I(end-1,:)];

%% calculate gradient angle
phi = atan2(Iy,Ix)*180/pi;

%% calculate length for horizontal / vertical edges
t = 8; %angle tolerance t
w_JNB = 3;
[row_idx, col_idx] = find(edges);
for k=1:length(row_idx)
    i = row_idx(k);
    j = col_idx(k);
    width_up=0; width_down=0;
    if (Ix(i,j) == 0 && Iy(i,j) == 0)
        continue; % not really neccesary
    end
    %% check for horizontal edge, gradient pointing upwards -> ~ 90°, ~ -270°
    if( abs(phi(i,j)+90) < t )
        min = 0;
        max = 0;
        for d = 1:m
            up = i-d;
            if (up < 1)
                width_up =  - 1;
                break;
            end
            if( I(up,j) <= I(up+1,j) ) %up+1 is max
              width_up = d - 1;
              max = I(up+1,j);
              break;
            end
        end
        for d = 1:m
            down = i+d;
            if (down > m)
                width_down =  - 1;
                break;
            end
            if( I(down,j) >= I(down-1,j) ) % down-1 is min
                width_down = d - 1;
                min = I(down-1,j);
                break;
            end
        end
        if(width_up ~= -1 && width_down ~= -1)
            widths_count = widths_count+1;
            phi2 = (phi(i,j)+90)*pi/180;
            edge_widths(i,j) = (width_up+width_down)/cos(phi2);
            slope = (max-min) / edge_widths(i,j);
            if (edge_widths(i,j) >= w_JNB)
                edge_widths(i,j) = edge_widths(i,j) - slope;
            end           
        end
    end  

    %% check for horizontal edge, gradient pointing downwards - -> ~ -90°, ~ 270°
    if( abs(phi(i,j)-90) < t ) 
        min = 0;
        max = 0;
        for d = 1:m
            up = i-d;
            if (up < 1)
                width_up =  - 1;
                break;
            end
            if( I(up,j) >= I(up+1,j) ) % up+1 is min
                width_up = d - 1;
                min = I(up+1,j);
                break;
            end
        end
        for d = 1:m
            down = i+d;
            if (down > m)
                width_down =  - 1;
                break;
            end
            if( I(down,j) <= I(down-1,j) ) %down-1 is max
                width_down = d - 1;
                max = I(down-1,j);
                break;
            end
        end
        if(width_up ~= -1 && width_down ~= -1)
            widths_count = widths_count+1;
            phi2 = (phi(i,j)-90) *pi/180;
            edge_widths(i,j) = (width_up+width_down)/cos(phi2);  
            slope = (max-min) / edge_widths(i,j);
            if (edge_widths(i,j) >= w_JNB)
                edge_widths(i,j) = edge_widths(i,j) - slope;
            end                     
        end
    end                     
end

%% calculate average edge widths in each block
avg_w = zeros(row_blocks,col_blocks);
for i=2:row_blocks-1 % skipping image borders 
    for j=2:col_blocks-1 
        block_row = (i-1)*BLOCKSIZE;
        block_col = (j-1)*BLOCKSIZE;
        block_widths = edge_widths(block_row+1:block_row+BLOCKSIZE,block_col+1:block_col+BLOCKSIZE);
        w_sum = sum(sum(block_widths));
        if ( w_sum >= THRESHOLD_W ) % enough widths found
            % calculate average width of the whole block
            avg_w(i,j) = w_sum / (sum(sum(block_widths ~= 0)));      
        end
    end
end
avg_w = avg_w(avg_w ~= 0);
avg_w = avg_w(:);
nr_of_used_blocks = ceil( numel(avg_w) * QUOTA_W );
if (nr_of_used_blocks == 0)
    s = 0;
    return;
end

avg_sorted = sort(avg_w);
sharpest_edges = avg_sorted(1:nr_of_used_blocks);

if (widths_count == 0)
    s=0;
else
    s = numel(sharpest_edges) / sum(sharpest_edges);
end

end



