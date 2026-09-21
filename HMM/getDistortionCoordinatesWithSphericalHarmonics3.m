function [r2PCS, D, rGrad] = getDistortionCoordinatesWithSphericalHarmonics3(rPCS, dicomHeader, imageSize)

    R0 = 0.25*1000;
    
    r0 = getSliceIsoCenter(dicomHeader);
    Q0 = convertVB17PatientPositionToGradientOrientation(dicomHeader.PatientPosition);
    
    T2 = [Q0 r0; 0 0 0 1];
    
    XYZ = T2\[rPCS; ones(1, size(rPCS, 2))];
    
    rGrad = XYZ(1:3,:);
    
    X = reshape(XYZ(1,:), imageSize);
    Y = reshape(XYZ(2,:), imageSize);
    Z = reshape(XYZ(3,:), imageSize);
    
    sphericalHarmonicsCoefficients = getSiemensLegendrePolynomialCoefficientsForSkyra();
    
    [ETA, D] = getDistortionField3(sphericalHarmonicsCoefficients, X, Y, Z, R0);
    
    X2 = X + ETA(:,:,:,1);
    Y2 = Y + ETA(:,:,:,2);
    Z2 = Z + ETA(:,:,:,3);
    
    XYZ2 = [X2(:)'; Y2(:)'; Z2(:)'; ones(1, prod(imageSize))];
    
    r2PCS = T2*XYZ2;
    
    r2PCS = r2PCS(1:3,:);
    
end