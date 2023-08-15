function [img_dose] = convertEBT3netODtoDose(EBT3, channel, calib, ...
    ctrl, bckg, ind)

counter = 0;
img_dose = [];

for i = ind % i = 1:length(EBT3) %  
    
    counter = counter + 1;
    
    % Convert pixel values into netOD, and then netOD into dose
    img_netOD = calculateNetOD(ctrl.(channel).PV_avg, ...
        EBT3{i}.(channel).img, bckg.(channel).PV_avg);
    temp_img_dose = model1(img_netOD, ...
        calib.(channel).mdl.Coefficients.Estimate(1), ...
        calib.(channel).mdl.Coefficients.Estimate(2), ...
        calib.(channel).n);
    img_dose = cat(3, img_dose, temp_img_dose);
    %     EBT3{i}.(channel).img_dose = img_dose;
    
    % Plot
    
%     figure()
%     hold on
%     imshow(img_netOD, [0 1], 'colormap', jet(4096)) % 5 Gy dose delivered
%     shading interp
%     c = colorbar;
%     c.Label.String = '\it{netOD}';
%     %     title(strrep(struct{i}.filename, '_', '-'))
%     %     leg = legend(['D_{max}=' num2str()]);
%     %     title(leg, struct{i}.filename)
%     set(gca, 'FontSize', 14)
% %     %     saveas(h, sprintf('%s.png', fn))


%     figure()
%     hold on
%     imshow(temp_img_dose, [0 6], 'colormap', jet(4096)) % 5 Gy dose delivered
%     shading interp
%     c = colorbar;
%     c.Label.String = 'Dose (Gy)';
% %     title(strrep(struct{i}.filename, '_', '-'))
% %     leg = legend(['D_{max}=' num2str()]);
% %     title(leg, struct{i}.filename)
%     set(gca, 'FontSize', 14)
% %     saveas(h, sprintf('%s.png', fn))
    
end

end