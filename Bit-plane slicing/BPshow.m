function fig = BPshow(BitPlanes)
figure()
for i = 1:length(BitPlanes)
    disp(strcat('Plotting Bit Plane',32,num2str(i)))
    subplot(2,4,i)
    imagesc(flipud(BitPlanes{i}));colormap('gray')
    title(strcat('Bit plane',32,num2str(i)))
end
disp('End of plotting. Wait for the figure to pop out')

end