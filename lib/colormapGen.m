%%
% generate color map
% input:
%  Cs            colors (color num x 3 matrix)
%  Is            color index pivot ([0, 1] color num vector)
%  csz           output colormap size
function cmap = colormapGen(Cs, Is, csz)
    cmap = zeros(csz,3);
    if any(Cs(:) > 1 & Cs(:) < 256), Cs = Cs / 255; end
    Is = Is * csz;
    Is(Is<=0) = 1;
    for i=1:length(Cs)-1
        step = Is(i+1)-Is(i);
        rgb = zeros(step+1,3);
        for j=1:3
            if (Cs(i+1,j)-Cs(i,j)) == 0
                rgb(:,j) = Cs(i,j);
            else
                rgb(:,j) = Cs(i,j):(Cs(i+1,j)-Cs(i,j))/step:Cs(i+1,j);
            end
        end
        cmap(Is(i):Is(i+1),:) = rgb;
    end
end
