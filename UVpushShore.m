function [U_final, V_final] = UVpushShore(Uin, Vin, dim)
%this function modifies a given velocity field (U,V) to free-slip conditions,
%so that shore points (i.e. ground points that have at least 1 ocean neighbor)
%acquire a velocity that is in a tangent direction to the shore at the
%point.

for i = 1:size(Uin,dim+1)
    
    if dim == 3
        %% for the full 3D (or 4D?) matrix:
        %pad with zeros:
        U = padarray(Uin(:,:,:,i),[1 1],0,'both');
        V = padarray(Vin(:,:,:,i),[1 1],0,'both');
        
        %ground indicator function:
        GroundInd = (U(:,:,:,1))==0; %gives 1 for ground, 0 for ocean, since U==0 only on the ground. Same borders for U and V so it doesn't matter which one is used for the calculation.
        
        
        
        %% detect shore points:
        %Detect the edges between land and ocean nodes by computing
        %the Laplacian with the 4 nearest neighbors:
        lapx_GInd_diff = diff(GroundInd(:,2:end,:),1,2) - diff(GroundInd(:,1:end-1,:),1,2);
        %if its value is negative, i is a shore point:
        lapx_GInd = lapx_GInd_diff<0; %lapx_GInd == 1 if it is a shore point in the x direction, 0 otherwise.
        %%%%%%%%%
        %Explanation:
        %for a given point i, Uind around this point in the x direction is
        %one of the following: (1,1,1), (1,1,0), (1,0,1), (1,0,0), (0,1,1),
        %(0,1,0), (0,0,1), (0,0,0).
        %The Laplacian calculates (Uind(i+1) - Uind(i)) - (Uind(i) - Uind(i-1)),
        %therefore will be negative when the i'th point is (1,1,0),
        %(0,1,1), (0,1,0) - i.e. when the i'th point is ground that has an
        %ocean neighbor in the x direction:
        %%%%%%%%%
        
        %same idea for y direction and diagonal directions:
        %y:
        lapy_GInd_diff = diff(GroundInd(2:end,:,:)) - diff(GroundInd(1:end-1,:,:));
        lapy_GInd = lapy_GInd_diff<0;
        %diag1:
        lapDIAG1_GInd_diff = GroundInd(1:end-2,1:end-2,:) + GroundInd(3:end, 3:end,:) - 2*GroundInd(2:end-1, 2:end-1,:);
        lapDIAG1_GInd = lapDIAG1_GInd_diff<0;
        %diag2:
        lapDIAG2_GInd_diff = GroundInd(3:end,1:end-2,:) + GroundInd(1:end-2, 3:end,:) - 2*GroundInd(2:end-1, 2:end-1,:);
        lapDIAG2_GInd = lapDIAG2_GInd_diff<0;
        
        %suffices that one of the laplacian indicator values is 1 so that the
        %point will be a shore point:
        Shore = sign(lapy_GInd(:,2:end-1,:) + lapx_GInd(2:end-1,:,:) + lapDIAG1_GInd + lapDIAG2_GInd);
        %Shore == 1 on shore points (ground points that have an ocean
        %border in some direction), and 0 on bulk ground and bulk ocean
        %points.
        
        %% Erase isolated shore/ocean points (with 3 out of 4 nearest neighbors ocean/shore respectively)
        %Identify isolated shore/ocean points:
        %isolated shore points:
        isolatedground_Ind = lapx_GInd_diff(2:end-1,:,:) + lapy_GInd_diff(:,2:end-1,:);
        IsolatedGround = (isolatedground_Ind == -3);
        
        %For isolated shore, change the velocity at these points to average
        %velocity of 3 neighboring ocean points:
        temp1 = U(1:end-2,:,:) + U(2:end-1,:,:) + U(3:end,:,:);
        U_IG = temp1(:,1:end-2,:) + temp1(:,2:end-1,:) + temp1(:,3:end,:);
        temp1 = V(1:end-2,:,:) + V(2:end-1,:,:) + V(3:end,:,:);
        V_IG = temp1(:,1:end-2,:) + temp1(:,2:end-1,:) + temp1(:,3:end,:);
        
        OceanInd = 1-GroundInd;
        temp = OceanInd(1:end-2,:,:) + OceanInd(2:end-1,:,:) + OceanInd(3:end,:,:);
        num_divide = (temp(:,1:end-2,:) + temp(:,2:end-1,:) + temp(:,3:end,:)) .* IsolatedGround;
        
        U_IG_switch = (U_IG.*IsolatedGround) ./ num_divide;
        V_IG_switch = (V_IG.*IsolatedGround) ./ num_divide;
        U_IG_switch(isnan(U_IG_switch)) = 0;
        V_IG_switch(isnan(V_IG_switch)) = 0;
        
        %and change to new shore without isolated points:
        Shore = Shore.*(1-IsolatedGround);
        GroundInd = GroundInd.*(1-padarray(IsolatedGround, [1,1],0,'both'));
        U = U + padarray(U_IG_switch, [1,1],0,'both');
        V = V + padarray(V_IG_switch, [1,1],0,'both');
        %For isolated ocean, change velocity to zero:
        
        %% Find inner corner points, that will later have a different rule for tangent:
        %inner corner points are ground points for which the sum of the 9
        %GroundInd values of them and their 8 nearest neighbors are 8:
        temp1 = GroundInd(1:end-2,:,:) + GroundInd(2:end-1,:,:) + GroundInd(3:end,:,:);
        temp2 = temp1(:,1:end-2,:) + temp1(:,2:end-1,:) + temp1(:,3:end,:);
        InnerCorner = (temp2 == 8) .* GroundInd(2:end-1, 2:end-1,:);
        
        
        %% calculate normal and tangent vectors to shore at each shore point:
        
        %For the shore points,
        %the derivative of GroundInd captures the component of the normal vector
        %in the direction of the derivative:
        %x direction:
        normalx_Uind = diff(GroundInd(:,1:end-1,:),1,2) + diff(GroundInd(:,2:end,:),1,2);
        normalx = normalx_Uind(2:end-1,:,:).*Shore;
        %y direction:
        normaly_Uind = diff(GroundInd(1:end-1,:,:)) + diff(GroundInd(2:end,:,:));
        normaly = normaly_Uind(:,2:end-1,:) .* Shore;
        
        %identify points that need a diagonal component:
        nonormalyet = (normalx==0 & normaly==0).*Shore;
        %add diagonal normal in diag1 direction:
        normalDIAG1_Uind = GroundInd(3:end, 3:end,:) - GroundInd(1:end-2,1:end-2,:);
        normalx2 = normalx + sqrt(2) * (normalDIAG1_Uind .* nonormalyet);
        normaly2 = normaly + sqrt(2) * (normalDIAG1_Uind .* nonormalyet);
        %add diagonal normal in diag2 direction:
        normalDIAG2_Uind = GroundInd(3:end, 1:end-2,:) - GroundInd(1:end-2,3:end,:);
        normalx3 = normalx2 - sqrt(2) * (normalDIAG2_Uind .* nonormalyet);
        normaly3 = normaly2 + sqrt(2) * (normalDIAG2_Uind .* nonormalyet);
        %%%%%%%
        %explanation:
        %i.e. in the x direction, a shore point has value GroundInd==1, with 
        %its neighbors either (1,1,0), (0,1,1), or (0,1,0). So the
        %derivative is 1/-1/0 depending on the relevant direction.
        %similar idea in other directions.
        %%%%%%%
        
        %normalize so that normal will be vector of length 1:
        %multiply by -1 (we actually calculated normal into ground instead of into ocean):
        Normalx = -normalx3 ./ sqrt(normalx3.^2 + normaly3.^2);
        Normaly = -normaly3 ./ sqrt(normalx3.^2 + normaly3.^2);
        
        %the tangent vector is just the 90° rotation of the normal vector:
        Tangentx = -Normaly;
        Tangenty = Normalx;
        
        %in inner corners, have the tangent in the same direction as the
        % normal to push particles out of the inner corners:
        Tangentx(InnerCorner==1) = Normalx(InnerCorner==1);
        Tangenty(InnerCorner==1) = Normaly(InnerCorner==1);
        
        %% identify relevant neighboring ocean velocity: follow normal to next grid point
        
        % check what neighbors need to be used for isolated shore points -
        % perhaps use nn in y direction for V and nn in x direction for U.
        
        %add comments to this code
        
        %take the normal vector
        neighbors_x = sign(-normalx3);
        %neighbors_x = sign(Normalx);
        neighbors_y = sign(-normaly3);
        %neighbors_y = sign(Normaly);
        
        %change size(Tangentx) to L or other variable...
        
        xindex = repmat((1:size(Tangentx,2)),size(Tangentx,1),1,size(Tangentx,3)) + neighbors_x;
        yindex = repmat((1:size(Tangentx,1))',1,size(Tangentx,2),size(Tangentx,3)) + neighbors_y;
        zindex = repmat(reshape(1:50,1,1,50), size(Tangentx,1), size(Tangentx,2), 1);
        
        xindex(xindex>size(Tangentx,2)) = size(Tangentx,2);
        yindex(yindex>size(Tangentx,1)) = size(Tangentx,1);
        xindex(xindex<1) = 1;
        yindex(yindex<1) = 1;
        
        xyzlist = [xindex(:), yindex(:), zindex(:)];
        ind = sub2ind(size(Tangentx), xyzlist(:,2), xyzlist(:,1), xyzlist(:,3));
        
        U1 = U(2:end-1,2:end-1,:,1);
        U1(isnan(U1)) = 0;
        U1_switch = reshape(U1(ind),size(U1)) .* Shore;
        
        V1 = V(2:end-1,2:end-1,:,1);
        V1(isnan(V1)) = 0;
        V1_switch = reshape(V1(ind),size(V1)) .* Shore;
        
        %length of velocity vector component that is parallel to tangent:
        %check if should not be U1_switch .* Tangentx + V1_switch .*
        %Tangenty; ???
        velLength = U1_switch .* Normaly + V1_switch .* Normalx;
        
        %switch velLength in inner corners to abs(velLength) so that the
        %velocity field will push particles out of the corners:
        %velLength(InnerCorner==1) = abs(velLength(InnerCorner==1));
        
        %x,y-components of new velocity:
        U1_new = abs(velLength) .* Normaly;
        V1_new = abs(velLength) .* Normalx;
        U1_new(isnan(U1_new)) = 0;
        V1_new(isnan(V1_new)) = 0;
        
        Vnew = V1 + V1_new;
        Unew = U1 + U1_new;
        
        changetonan = sign(Vnew).^2 + sign(Unew).^2;
        Unew(changetonan == 0) = nan;
        Vnew(changetonan == 0) = nan;
        
        
        %% sanity check:
        % figure; imagesc(Uind(2:end-1, 2:end-1,1))
        % hold on; quiver(-normalx3(:,:,1), -normaly3(:,:,1),'m')
        % caxis([0 10])
        % hold on; quiver(Tangentx(:,:,1), Tangenty(:,:,1),'k','LineWidth', 2)
        % hold on; quiver(V1(:,:,1), U1(:,:,1),'r','LineWidth', 2)
        % hold on; quiver(V1_switch(:,:,1), U1_switch(:,:,1),'g','LineWidth', 2)
        % hold on; quiver(V1_final(:,:,1), U1_final(:,:,1),'w','LineWidth', 1)
        
        %%
        %opposite??
        U_final(:,:,:,i) = Unew;
        V_final(:,:,:,i) = Vnew;
        
        
    elseif dim == 2
        %% for a 2D flow:
        %pad with zeros:
        U = padarray(Uin(:,:,i),[1 1],0,'both');
        V = padarray(Vin(:,:,i),[1 1],0,'both');
        
        Uindsmall = isnan(U);
        
        lapx_Uindsmall = diff(Uindsmall(:,2:end),1,2) - diff(Uindsmall(:,1:end-1),1,2);
        lapx_Uindsmall = lapx_Uindsmall<0;
        
        lapy_Uindsmall = diff(Uindsmall(2:end,:)) - diff(Uindsmall(1:end-1,:));
        lapy_Uindsmall = lapy_Uindsmall<0;
        
        lapDIAG1_Uindsmall = Uindsmall(1:end-2,1:end-2) + Uindsmall(3:end, 3:end) - 2*Uindsmall(2:end-1, 2:end-1);
        lapDIAG1_Uindsmall = lapDIAG1_Uindsmall<0;
        
        lapDIAG2_Uindsmall = Uindsmall(3:end,1:end-2) + Uindsmall(1:end-2, 3:end) - 2*Uindsmall(2:end-1, 2:end-1);
        lapDIAG2_Uindsmall = lapDIAG2_Uindsmall<0;
        
        Shore = lapy_Uindsmall(:,2:end-1) + lapx_Uindsmall(2:end-1,:) + lapDIAG1_Uindsmall + lapDIAG2_Uindsmall;
        Shore = sign(Shore);
        
        
        
        normalx_Uindsmall = diff(Uindsmall(:,1:end-1),1,2) + diff(Uindsmall(:,2:end),1,2);
        normalx = normalx_Uindsmall(2:end-1,:).*Shore;
        
        normaly_Uindsmall = diff(Uindsmall(1:end-1,:)) + diff(Uindsmall(2:end,:));
        normaly = normaly_Uindsmall(:,2:end-1) .* Shore;
        
        nonormalyet = (normalx==0 & normaly==0).*Shore;
        
        normalDIAG1_Uindsmall = Uindsmall(3:end, 3:end) - Uindsmall(1:end-2,1:end-2);
        normalx2 = normalx + sqrt(2) * (normalDIAG1_Uindsmall .* nonormalyet);
        normaly2 = normaly + sqrt(2) * (normalDIAG1_Uindsmall .* nonormalyet);
        
        normalDIAG2_Uindsmall = Uindsmall(3:end, 1:end-2) - Uindsmall(1:end-2,3:end);
        normalx3 = normalx2 - sqrt(2) * (normalDIAG2_Uindsmall .* nonormalyet);
        normaly3 = normaly2 + sqrt(2) * (normalDIAG2_Uindsmall .* nonormalyet);
        
        % figure; imagesc(Uindsmall(2:end-1, 2:end-1))
        % hold on; quiver(-normalx3, -normaly3,'m','LineWidth', 6)
        
        Normalx = -normalx3 ./ sqrt(normalx3.^2 + normaly3.^2);
        Normaly = -normaly3 ./ sqrt(normalx3.^2 + normaly3.^2);
        
        Tangentx = -Normaly./sqrt(Normalx.^2 + Normaly.^2);
        Tangenty = Normalx./sqrt(Normalx.^2 + Normaly.^2);
        Tangentx(isnan(Tangentx)) = 0;
        Tangenty(isnan(Tangenty)) = 0;
        
        %% identify relevant neighboring ocean velocity: follow normal to next grid point
        
        neighbors_x = sign(-normalx3);
        neighbors_y = sign(-normaly3);
        
        xindex = repmat((1:19),19,1) + neighbors_x;
        yindex = repmat((1:19)',1,19) + neighbors_y;
        
        xindex(xindex>19) = 19;
        yindex(yindex>19) = 19;
        xindex(xindex<1) = 1;
        yindex(yindex<1) = 1;
        
        xylist = [xindex(:), yindex(:)];
        ind = sub2ind(size(Tangentx), xylist(:,2), xylist(:,1));
        
        Usmall = U(101:119, 101:119, 1);
        Usmall(isnan(Usmall)) = 0;
        Usmall_switch = reshape(Usmall(ind),size(Usmall)) .* Shore;
        
        
        Vsmall = V(101:119, 101:119, 1);
        Vsmall(isnan(Vsmall)) = 0;
        Vsmall_switch = reshape(Vsmall(ind),size(Vsmall))  .* Shore;
        
        %length of velocity vector component that is parallel to tangent:
        velLength = Usmall_switch .* Tangenty + Vsmall_switch .* Tangentx;
        
        %x,y-components of new velocity:
        Usmall_new = velLength .* Tangenty;
        Vsmall_new = velLength .* Tangentx;
        
        Vnew = Vsmall + Vsmall_new;
        Unew = Usmall + Usmall_new;
        
        changetonan = sign(Vnew + Unew).^2;
        Unew(changetonan == 0) = nan;
        Vnew(changetonan == 0) = nan;
        
        
        % %sanity check:
        % figure; imagesc(Uindsmall(2:end-1, 2:end-1,1))
        % caxis([0 10])
        % hold on; quiver(-normalx3(:,:), -normaly3(:,:),'m','LineWidth', 2)
        % hold on; quiver(Tangentx, Tangenty,'k','LineWidth', 1)
        % hold on; quiver(Vsmall, Usmall,'r','LineWidth', 2)
        % hold on; quiver(Vsmall_switch, Usmall_switch,'g','LineWidth', 2)
        % hold on; quiver(Vsmall_final, Usmall_final,'w','LineWidth', 3)
        
        U_final(:,:,i) = Unew;
        V_final(:,:,i) = Vnew;
    end
end