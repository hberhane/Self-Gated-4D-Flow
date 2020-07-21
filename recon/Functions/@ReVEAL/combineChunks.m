function combineChunks(obj, dim)

switch dim
    
    % Combine along frequency encode (RO) dimension
    case 1
        
        N(1) = obj.data.chunkInds(end,2);
        N(2) = size(obj.data.chunkOutputs(1).xHatb,2);
        N(3) = size(obj.data.chunkOutputs(1).xHatb,3);
        N(4) = size(obj.data.chunkOutputs(1).xHatb,4);
        
        xHatb = zeros(N);
        xHatx = zeros(N);
        xHaty = zeros(N);
        xHatz = zeros(N);
        
        thetaX = zeros(N);
        thetaY = zeros(N);
        thetaZ = zeros(N);
        
        vX = zeros(N);
        vY = zeros(N);
        vZ = zeros(N);
        
        xHatb(obj.data.chunkInds(1,1):obj.data.chunkInds(1,2),:,:,:) = obj.data.chunkOutputs(1).xHatb;
        xHatx(obj.data.chunkInds(1,1):obj.data.chunkInds(1,2),:,:,:) = obj.data.chunkOutputs(1).xHatx;
        xHaty(obj.data.chunkInds(1,1):obj.data.chunkInds(1,2),:,:,:) = obj.data.chunkOutputs(1).xHaty;
        xHatz(obj.data.chunkInds(1,1):obj.data.chunkInds(1,2),:,:,:) = obj.data.chunkOutputs(1).xHatz;

        thetaX(obj.data.chunkInds(1,1):obj.data.chunkInds(1,2),:,:,:) = obj.data.chunkOutputs(1).thetaX;
        thetaY(obj.data.chunkInds(1,1):obj.data.chunkInds(1,2),:,:,:) = obj.data.chunkOutputs(1).thetaY;
        thetaZ(obj.data.chunkInds(1,1):obj.data.chunkInds(1,2),:,:,:) = obj.data.chunkOutputs(1).thetaZ;
        
        vX(obj.data.chunkInds(1,1):obj.data.chunkInds(1,2),:,:,:) = obj.data.chunkOutputs(1).vX;
        vY(obj.data.chunkInds(1,1):obj.data.chunkInds(1,2),:,:,:) = obj.data.chunkOutputs(1).vY;
        vZ(obj.data.chunkInds(1,1):obj.data.chunkInds(1,2),:,:,:) = obj.data.chunkOutputs(1).vZ;        

        for i = 2 : length(obj.data.chunkOutputs)
            
            ovLap = obj.data.chunkInds(i-1,2) - obj.data.chunkInds(i,1) + 1;
            
            xHatb(obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2),:,:,:) = 0.5 * (obj.data.chunkOutputs(i-1).xHatb((end-ovLap+1):end,:,:,:) + obj.data.chunkOutputs(i).xHatb(1:ovLap,:,:,:));
            xHatb(obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2),:,:,:) = obj.data.chunkOutputs(i).xHatb(ovLap+1:end,:,:,:);

            xHatx(obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2),:,:,:) = 0.5 * (obj.data.chunkOutputs(i-1).xHatx((end-ovLap+1):end,:,:,:) + obj.data.chunkOutputs(i).xHatx(1:ovLap,:,:,:));
            xHatx(obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2),:,:,:) = obj.data.chunkOutputs(i).xHatx(ovLap+1:end,:,:,:);
            
            xHaty(obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2),:,:,:) = 0.5 * (obj.data.chunkOutputs(i-1).xHaty((end-ovLap+1):end,:,:,:) + obj.data.chunkOutputs(i).xHaty(1:ovLap,:,:,:));
            xHaty(obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2),:,:,:) = obj.data.chunkOutputs(i).xHaty(ovLap+1:end,:,:,:);            
            
            xHatz(obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2),:,:,:) = 0.5 * (obj.data.chunkOutputs(i-1).xHatz((end-ovLap+1):end,:,:,:) + obj.data.chunkOutputs(i).xHatz(1:ovLap,:,:,:));
            xHatz(obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2),:,:,:) = obj.data.chunkOutputs(i).xHatz(ovLap+1:end,:,:,:);
            
            

            thetaX(obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2),:,:,:) = 0.5 * (obj.data.chunkOutputs(i-1).thetaX((end-ovLap+1):end,:,:,:) + obj.data.chunkOutputs(i).thetaX(1:ovLap,:,:,:));
            thetaX(obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2),:,:,:) = obj.data.chunkOutputs(i).thetaX(ovLap+1:end,:,:,:);
            
            thetaY(obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2),:,:,:) = 0.5 * (obj.data.chunkOutputs(i-1).thetaY((end-ovLap+1):end,:,:,:) + obj.data.chunkOutputs(i).thetaY(1:ovLap,:,:,:));
            thetaY(obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2),:,:,:) = obj.data.chunkOutputs(i).thetaY(ovLap+1:end,:,:,:);          
            
            thetaZ(obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2),:,:,:) = 0.5 * (obj.data.chunkOutputs(i-1).thetaZ((end-ovLap+1):end,:,:,:) + obj.data.chunkOutputs(i).thetaZ(1:ovLap,:,:,:));
            thetaZ(obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2),:,:,:) = obj.data.chunkOutputs(i).thetaZ(ovLap+1:end,:,:,:);
            
            
            vX(obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2),:,:,:) = 0.5 * (obj.data.chunkOutputs(i-1).vX((end-ovLap+1):end,:,:,:) + obj.data.chunkOutputs(i).vX(1:ovLap,:,:,:));
            vX(obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2),:,:,:) = obj.data.chunkOutputs(i).vX(ovLap+1:end,:,:,:);
            
            vY(obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2),:,:,:) = 0.5 * (obj.data.chunkOutputs(i-1).vY((end-ovLap+1):end,:,:,:) + obj.data.chunkOutputs(i).vY(1:ovLap,:,:,:));
            vY(obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2),:,:,:) = obj.data.chunkOutputs(i).vY(ovLap+1:end,:,:,:);          
            
            vZ(obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2),:,:,:) = 0.5 * (obj.data.chunkOutputs(i-1).vZ((end-ovLap+1):end,:,:,:) + obj.data.chunkOutputs(i).vZ(1:ovLap,:,:,:));
            vZ(obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2),:,:,:) = obj.data.chunkOutputs(i).vZ(ovLap+1:end,:,:,:);
            
        end
        
        obj.outputs.xHatb = xHatb;
        obj.outputs.xHatx = xHatx;
        obj.outputs.xHaty = xHaty;
        obj.outputs.xHatz = xHatz;
        
        obj.outputs.thetaX = thetaX;
        obj.outputs.thetaY = thetaY;
        obj.outputs.thetaZ = thetaZ;
        
        obj.outputs.vX = vX;
        obj.outputs.vY = vY;
        obj.outputs.vZ = vZ;
    
    % Combine along cardiac dimension
    case 5
        
        N(1) = size(obj.data.chunkOutputs(1).xHatb,1);
        N(2) = size(obj.data.chunkOutputs(1).xHatb,2);
        N(3) = size(obj.data.chunkOutputs(1).xHatb,3);
        N(4) = obj.data.chunkInds(end,2);

        xHatb = zeros(N);
        xHatx = zeros(N);
        xHaty = zeros(N);
        xHatz = zeros(N);
        
        thetaX = zeros(N);
        thetaY = zeros(N);
        thetaZ = zeros(N);
        
        vX = zeros(N);
        vY = zeros(N);
        vZ = zeros(N);
        
        xHatb(:,:,:,obj.data.chunkInds(1,1):obj.data.chunkInds(1,2)) = obj.data.chunkOutputs(1).xHatb;
        xHatx(:,:,:,obj.data.chunkInds(1,1):obj.data.chunkInds(1,2)) = obj.data.chunkOutputs(1).xHatx;
        xHaty(:,:,:,obj.data.chunkInds(1,1):obj.data.chunkInds(1,2)) = obj.data.chunkOutputs(1).xHaty;
        xHatz(:,:,:,obj.data.chunkInds(1,1):obj.data.chunkInds(1,2)) = obj.data.chunkOutputs(1).xHatz;
        
        
        
        for i = 2 : length(obj.data.chunkOutputs)
            
            ovLap = obj.data.chunkInds(i-1,2) - obj.data.chunkInds(i,1) + 1;
            
            xHatb(:,:,:,obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2)) = 0.5 * (obj.data.chunkOutputs(i-1).xHatb(:,:,:,(end-ovLap+1):end) + obj.data.chunkOutputs(i).xHatb(:,:,:,1:ovLap));
            xHatb(:,:,:,obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2)) = obj.data.chunkOutputs(i).xHatb(:,:,:,ovLap+1:end);

            xHatx(:,:,:,obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2)) = 0.5 * (obj.data.chunkOutputs(i-1).xHatx(:,:,:,(end-ovLap+1):end) + obj.data.chunkOutputs(i).xHatx(:,:,:,1:ovLap));
            xHatx(:,:,:,obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2)) = obj.data.chunkOutputs(i).xHatx(:,:,:,ovLap+1:end);
            
            xHaty(:,:,:,obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2)) = 0.5 * (obj.data.chunkOutputs(i-1).xHaty(:,:,:,(end-ovLap+1):end) + obj.data.chunkOutputs(i).xHaty(:,:,:,1:ovLap));
            xHaty(:,:,:,obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2)) = obj.data.chunkOutputs(i).xHaty(:,:,:,ovLap+1:end);            
            
            xHatz(:,:,:,obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2)) = 0.5 * (obj.data.chunkOutputs(i-1).xHatz(:,:,:,(end-ovLap+1):end) + obj.data.chunkOutputs(i).xHatz(:,:,:,1:ovLap));
            xHatz(:,:,:,obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2)) = obj.data.chunkOutputs(i).xHatz(:,:,:,ovLap+1:end);
            
            
            thetaX(:,:,:,obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2)) = 0.5 * (obj.data.chunkOutputs(i-1).thetaX(:,:,:,(end-ovLap+1):end) + obj.data.chunkOutputs(i).thetaX(:,:,:,1:ovLap));
            thetaX(:,:,:,obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2)) = obj.data.chunkOutputs(i).thetaX(:,:,:,ovLap+1:end);            

            thetaY(:,:,:,obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2)) = 0.5 * (obj.data.chunkOutputs(i-1).thetaY(:,:,:,(end-ovLap+1):end) + obj.data.chunkOutputs(i).thetaY(:,:,:,1:ovLap));
            thetaY(:,:,:,obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2)) = obj.data.chunkOutputs(i).thetaY(:,:,:,ovLap+1:end);
            
            thetaZ(:,:,:,obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2)) = 0.5 * (obj.data.chunkOutputs(i-1).thetaZ(:,:,:,(end-ovLap+1):end) + obj.data.chunkOutputs(i).thetaZ(:,:,:,1:ovLap));
            thetaZ(:,:,:,obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2)) = obj.data.chunkOutputs(i).thetaZ(:,:,:,ovLap+1:end);
            
            
            vX(:,:,:,obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2)) = 0.5 * (obj.data.chunkOutputs(i-1).vX(:,:,:,(end-ovLap+1):end) + obj.data.chunkOutputs(i).vX(:,:,:,1:ovLap));
            vX(:,:,:,obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2)) = obj.data.chunkOutputs(i).vX(:,:,:,ovLap+1:end);
            
            vY(:,:,:,obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2)) = 0.5 * (obj.data.chunkOutputs(i-1).vY(:,:,:,(end-ovLap+1):end) + obj.data.chunkOutputs(i).vY(:,:,:,1:ovLap));
            vY(:,:,:,obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2)) = obj.data.chunkOutputs(i).vY(:,:,:,ovLap+1:end);        
            
            vZ(:,:,:,obj.data.chunkInds(i,1):obj.data.chunkInds(i-1,2)) = 0.5 * (obj.data.chunkOutputs(i-1).vZ(:,:,:,(end-ovLap+1):end) + obj.data.chunkOutputs(i).vZ(:,:,:,1:ovLap));
            vZ(:,:,:,obj.data.chunkInds(i-1,2)+1:obj.data.chunkInds(i,2)) = obj.data.chunkOutputs(i).vZ(:,:,:,ovLap+1:end);            
        end

obj.outputs.xHatb = xHatb;
obj.outputs.xHatx = xHatx;
obj.outputs.xHaty = xHaty;
obj.outputs.xHatz = xHatz;

obj.outputs.thetaX = thetaX;
obj.outputs.thetaY = thetaY;
obj.outputs.thetaZ = thetaZ;

obj.outputs.vX = vX;
obj.outputs.vY = vY;
obj.outputs.vZ = vZ;

end

