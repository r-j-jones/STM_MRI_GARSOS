function [ETA, D, OMEGA] = getDistortionField3(coefficentTable, X, Y, Z, R0)

R = sqrt(X.^2 + Y.^2 + Z.^2);
OMEGA = R <= R0;
THETA = acos(Z./R);
THETA(R == 0) = 0;
PHI = atan2(Y, X);
COSTHETA = cos(THETA);

ETA = zeros([size(X) 3]);

degrees = sort(unique(cell2mat(coefficentTable(:,2))),'ascend');
nDegrees = numel(degrees);

pg = Progressor('Generating distortion field.');

cTerm = 0;

for iDegree = 1:nDegrees
    
    l = degrees(iDegree);
    degreeCoefficentTable = coefficentTable(cell2mat(coefficentTable(:,2)) == l, :);

    Rn = (R/R0).^l;
    
    orders = sort(unique(cell2mat(degreeCoefficentTable(:,3))),'ascend');
    nOrders = numel(orders);
    
    for iOrder = 1:nOrders
        
        m = orders(iOrder);        
        orderCoefficentTable = degreeCoefficentTable(cell2mat(degreeCoefficentTable(:,3)) == m, :);
        
        if m == 0            
            P = Rn.*legendreP_mod20130508(l, m, COSTHETA);            
        else        
            P = Rn.*((-1)^m)*sqrt(((l + 0.5)*factorial(l - m))/factorial(l + m)).*legendreP_mod20130508(l, m, COSTHETA);        
        end
        
        terms = sort(unique(cell2mat(orderCoefficentTable(:,1))),'ascend');
        nTerms = numel(terms);
        
        for iTerm = 1:nTerms
            
            t = terms(iTerm);
            termCoefficentTable = orderCoefficentTable(cell2mat(orderCoefficentTable(:,1)) == t, :);
            
            switch t
                case 'A'
                    T = P.*cos(m*PHI);
                case 'B'
                    T = P.*sin(m*PHI);
            end
            
            dimensions = sort(unique(cell2mat(termCoefficentTable(:,5))),'ascend');
            nDimensions = numel(dimensions);
            
            for iDimension = 1:nDimensions
                d = dimensions(iDimension);
                
                C = T.*termCoefficentTable{cell2mat(termCoefficentTable(:,5)) == d, 4};
                
                di = double(d) - double('x') + 1;
                
                ETA(:, :, :, di) = ETA(:, :, :, di) + C;
                
%                 fprintf('%s(%s,%.0f,%.0f)\n',t,d,l,m);
                
                cTerm = cTerm + 1;
                pg.setProgress(cTerm/size(coefficentTable,1));
            end
            
        end
        
    end
end

ETA(:,:,:,1) = ETA(:,:,:,1)*R0;
ETA(:,:,:,2) = ETA(:,:,:,2)*R0;
ETA(:,:,:,3) = ETA(:,:,:,3)*R0;

dXdx = diff(ETA(:,:,:,1),1,1)./diff(X,1,1);
dYdy = diff(ETA(:,:,:,2),1,2)./diff(Y,1,2);
dZdz = diff(ETA(:,:,:,3),1,3)./diff(Z,1,3);

dXdx = (cat(1,dXdx(1,:,:),dXdx) + cat(1,dXdx,dXdx(end,:,:)))/2;
dYdy = (cat(2,dYdy(:,1,:),dYdy) + cat(2,dYdy,dYdy(:,end,:)))/2;
dZdz = (cat(3,dZdz(:,:,1),dZdz) + cat(3,dZdz,dZdz(:,:,end)))/2;

D = (1 - dXdx).*(1 - dYdy).*(1 - dZdz);


