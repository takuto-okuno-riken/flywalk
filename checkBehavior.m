% checking registered calcium imaging histogram, etc.

function checkBehavior
    frames = 3384; % frame num of calcium imaging
    offset = 104.654892; % ms, start time of calcium imaging
    TR = 1000 / 1.879714; % ms, retrieval time of calcium imaging

    % check movement speed
    listing = dir('behavior/sub-fly-*.mat');
    MS = [];
    for i=1:length(listing)
        folder = listing(i).folder;
        name = listing(i).name;

        load([folder '/' name]);
        MS = [MS; MovementSpeed]; % 50 Hz, 20ms

        figure; plot(Position(1,:), Position(2,:)); title(['trajectory of ' name]);
    end
    spm = mean(MS,2);
    spsd = std(MS,1,2);
    figure; plot(spm); ylim([0 0.03]); title('mean speed in each fly');
    figure; plot(spsd); ylim([0 0.03]); title('std speed in each fly');

    % get mean speed per TR
    MScal = zeros(size(MS,1),frames,'single');
    for i=1:frames
        st = ceil((offset + TR*(i-1)) / 20);
        ed = floor((offset + TR*i) / 20);
        if ed > size(MS,2), ed = size(MS,2); end
        MScal(:,i) = mean(MS(:,st:ed),2);
    end
    P = prctile(MS,10,2); % find bottom 10% speed
    bsp = mean(P);
    th = mean(spsd);
    MBcal = MScal;
    MBcal(MBcal>th) = th; MBcal = MBcal - bsp; MBcal(MBcal<0) = 0;
    MBcal = MBcal / (th - bsp); % normalize [0, 1]

    figure; plot(MScal'); title('speed per TR in each fly');
    figure; plot(MBcal'); title('moving behavior per TR in each fly');

    save('data/behavior.mat','MScal','MBcal');
end
