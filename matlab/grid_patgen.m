% pushbroom GRID  pattern generator
% dmd dimensions: 1140 x 912 
% grid spacing: 16px
% off pixels (0) correspond to the DMD mirror reflecting light into the
% camera, on pixels (255) mean that mirror has no contribution to the image

rows = 1140;    % rows in DMD
cols = 912;     % columns in DMD
tSize = 14;      % grid tile size in pixels

gRows = floor(rows / tSize);   % number of tiles vertically
gCols = floor(cols / tSize);   % number of tiles horizontally

rowActual = gRows*tSize;
colActual = gCols*tSize;

% now build pyramid + all black and all white patterns
% and prepend them to the shifting grid patterns
[P,M,N]=buildPatterns(6,0);
    
imageCount = 0; % for file naming

for imageCount=1:size(P,3)
    % write to disk
    I=P(:,:,imageCount);
    I=repmat(I,[1,1,3]);
    if imageCount<10
        imwrite(I, ['pyramidpat/pattern00' num2str(imageCount) '.bmp']);
    elseif imageCount>=10 && imageCount<100
        imwrite(I, ['pyramidpat/pattern0' num2str(imageCount) '.bmp']);
    else
        imwrite(I, ['pyramidpat/pattern' num2str(imageCount) '.bmp']);
    end
end

 for tr=1:tSize
    for tc=1:tSize
        pat = repmat(zeros(tSize,tSize),1); % generate the 16x16 tile
        pat(tr,tc)=255;
        I=repmat(pat,gRows,gCols); % create full size image
        
        I=[I;zeros(rows-rowActual,colActual)]; % pad with 4 rows at bottom
        if colActual ~= cols
            I=[I zeros(rows,cols-colActual)]; % pad with 4 rows at right
        end
      
        I=repmat(I,[1,1,3]);
        
        % write to disk
        if imageCount<10
            imwrite(I, ['pyramidpat/pattern00' num2str(imageCount) '.bmp']);
        elseif imageCount>=10 && imageCount<100
            imwrite(I, ['pyramidpat/pattern0' num2str(imageCount) '.bmp']);
        else
            imwrite(I, ['pyramidpat/pattern' num2str(imageCount) '.bmp']);
        end
        
        imageCount=imageCount+1
   end
end

