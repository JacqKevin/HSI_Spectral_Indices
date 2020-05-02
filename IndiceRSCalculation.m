function [C,label] = IndiceRSCalculation(M,RGB,d,wl,indx)
% Fonction qui permet de calculer les indices de la littérature en remote
% sensing.
% INPUT :
%           M : cube hyperspectral
%           wl : longueur d'onde
% OUTPUT :
%           C : carte des indices
%           label : indices calculés

if mean(RGB(:))>1000
    RGB=RGB/10000;
end
RGB=RGB*(0.5/mean(RGB(:)));

map_red=cat(2,[0.5:0.001:1;zeros(1,501);zeros(1,501)],[ones(1,1001);0:0.001:1;0:0.001:1])';
map_green=cat(2,[zeros(1,501);0.5:0.001:1;zeros(1,501)],[0:0.001:1;ones(1,1001);0:0.001:1])';
map_orange=[ones(1,1001);1:-0.0005:.5;zeros(1,1001)]';
map_blue=cat(2,[zeros(1,501);zeros(1,501);0.5:0.001:1],[0:0.001:1;0:0.001:1;ones(1,1001)])';
p_r=[1 0 0];
p_g=[0 1 0];
p_b=[0 0 1];
p_o=[1 1 0];

M=M/10000;

%% Chlorophyll
% Check the wavelengths range
iter=0;
if mean(wl)<1000
    
    index_vnir={'Normalized Difference Vegetation Index (NDVI)',...
        'Normalized difference water index (NDWI)',...
        'Modified Soil-adjusted Vegetation Index (MSAVI2)',...
        'Modified Triangular Vegetation Index (MTVI2)',...
        'Visible Atmospherically Resistant Index (VARI)',...
        'triangular greenness index (TGI)',...
        'Enhanced Vegetation Index (EVI)',...
        'Soil-adjusted Vegetation Index (SAVI)',...
        'Global Environmental Monitoring Index (GEMI)',...
        'Green Normalized Difference Vegetation Index (GNDVI)',...
        'Iron Oxide (IO)','SR680 Chlorophylls','SR705 Chlorophylls',...
        'ND680 Chlorophylls','ND705 Chlorophylls','mSR705 Chlorophylls',...
        'mND705 Chlorophylls','PRI Chlorophylls/Carotenoids',...
        'SIPI Chlorophylls/Carotenoids','PSRI Chlorophylls/Carotenoids',...
        'R/G Anthocyanin'};
    associated_vnir_map={map_green;map_blue;map_green;map_green;...
        map_green;map_green;map_green;map_green;map_green;map_green;...
        map_red;map_green;map_green;map_green;map_green;map_green;...
        map_green;map_green;map_green;map_green;map_orange};
    associated_vnir_p={p_g;p_b;p_g;p_g;p_g;p_g;p_g;p_g;p_g;p_g;p_r;...
        p_g;p_g;p_g;p_g;p_g;p_g;p_g;p_g;p_g;p_o};
    
    if nargin<5
        [indx,~] = listdlg('ListString',index_vnir);
    end
    
    [~,Blue]=find(abs(wl-470)==min(abs(wl-470)));
    [~,Green]=find(abs(wl-550)==min(abs(wl-550)));
    [~,Red]=find(abs(wl-650)==min(abs(wl-650)));
    [~,NIR]=find(abs(wl-860)==min(abs(wl-860)));
    
    if sum(indx==1) % NDVI Normalized Difference Vegetation Index
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Red(1))));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==2) % Normalized difference water index NDWI
        iter=iter+1;
        
        IM=(squeeze(M(:,:,Green(1)))-squeeze(M(:,:,NIR(1))))./...
            (squeeze(M(:,:,Green(1)))+squeeze(M(:,:,NIR(1))));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==3) % MSAVI2 Modified Soil-adjusted Vegetation Index
        iter=iter+1;
        
        IM=(2*squeeze(M(:,:,NIR(1)))+1-...
            sqrt((2*squeeze(M(:,:,NIR(1)))+1).^2-...
            8*(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))))/2;
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==4) % MTVI2 Modified Triangular Vegetation Index
        iter=iter+1;
        
        IM=(1.5*(1.2*(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Green(1))))-...
            2.5*(squeeze(M(:,:,Red(1)))-squeeze(M(:,:,Green(1))))))./...
            sqrt((2*squeeze(M(:,:,NIR(1)))+1).^2-...
            (6*squeeze(M(:,:,NIR(1)))-5*sqrt(squeeze(M(:,:,Red(1)))))-0.5);
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==5) % VARI Visible Atmospherically Resistant Index
        iter=iter+1;
        
        IM=(squeeze(M(:,:,Green(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,Green(1)))+squeeze(M(:,:,Red(1)))-...
            squeeze(M(:,:,Blue(1))));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==6) % TGI triangular greenness index
        iter=iter+1;
        
        IM=(120*(squeeze(M(:,:,Red(1)))-squeeze(M(:,:,Blue(1))))-...
            190*(squeeze(M(:,:,Red(1)))-squeeze(M(:,:,Green(1)))))/2;
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==7) % EVI  Enhanced Vegetation Index
        iter=iter+1;
        
        IM=2.5*(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1)))+6*squeeze(M(:,:,Red(1)))-...
            7.5*squeeze(M(:,:,Blue(1)))+1);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==8) % SAVI Soil-adjusted Vegetation Index
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Red(1)))+1)*(1+1);
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==9) % Global Environmental Monitoring Index (GEMI)
        iter=iter+1;
        
        eta=(2*(squeeze(M(:,:,NIR(1))).^2-squeeze(M(:,:,Red(1))).^2)+...
            1.5*squeeze(M(:,:,NIR(1)))+0.5*squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Red(1)))+0.5);
        IM=eta.*(1-0.25*eta)-((squeeze(M(:,:,Red(1)))-0.125)./...
            (1-squeeze(M(:,:,Red(1)))));
        clear eta
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==10) % Green Normalized Difference Vegetation Index (GNDVI)
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Green(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Green(1))));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==11) % Iron Oxide (IO)
        iter=iter+1;
        
        IM=squeeze(M(:,:,Red(1)))./squeeze(M(:,:,Blue(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==12) % SR680 Chlorophylls
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[800 680];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==13) % SR705 Chlorophylls
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[750 705];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==14) % ND680 Chlorophylls
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[800 680];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==15) % ND705 Chlorophylls
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[750 705];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==16) % mSR705 Chlorophylls
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[750 445 705];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind3(1)))-squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==17) % mND705 Chlorophylls
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[750 705 445];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1)))-...
            2*squeeze(M(:,:,ind3(1))));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==18) % PRI Chlorophylls/Carotenoids
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[531 570];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==19) % SIPI Chlorophylls/Carotenoids
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[800 445 680];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind3(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==20) % PSRI Chlorophylls/Carotenoids
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[680 500 750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind3(1))));
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(indx==21) % R/G Anthocyanin
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[699 600 500];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=sum(squeeze(M(:,:,ind2(1):ind1(1))),3)./...
            sum(squeeze(M(:,:,ind3(1):ind2(1))),3);
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    for i=1:size(C,3)
        ha=[];
        IM=squeeze(C(:,:,i));
        figure;
        ha(1)=subplot(3,1,1);
        imagesc(d,d(1:size(C,1)),IM)
        caxis([nanmean(IM(:))-3*nanstd(IM(:)) nanmean(IM(:))+3*nanstd(IM(:))])
        colorbar
        map=associated_vnir_map{indx(i)};
        colormap(map)
        xlabel('Depth (cm)')
        ylabel('Width (cm)')
        set(gca,'fontsize',14)
        title(index_vnir{indx(i)})
        ha(2)=subplot(3,1,2);
        imagesc(d,d(1:size(IM,1)),RGB)
        colorbar
        xlabel('Depth (cm)')
        ylabel('Width (cm)')
        set(gca,'fontsize',14)
        ha(3)=subplot(313);
        plot(d,nanmean(IM(round(0.4*size(IM,1):0.6*size(IM,1)),:),1),'color',associated_vnir_p{indx(i)})
        grid on
        xlabel('Depth (cm)')
        ylabel('SI')
        set(gca,'fontsize',14)
        colorbar
        xlim([min(d) max(d)])
        linkaxes(ha,'x')
    end
    
    figure;
    imagesc(corr(reshape(C,[],size(C,3))))
    colorbar
    colormap(jet)
    xtickangle(45)
    set(gca,'fontsize',14,'xtick',1:length(label),'xticklabel',label,...
        'ytick',1:length(label),'yticklabel',label)
end

if mean(wl)>1000
    
    index_swir={'CM','Sand moisture index (SMI)',...
'SWIR Fine particles Index (FI)'};
    if nargin<5
        [indx,~] = listdlg('ListString',index_swir);
    end
    
    [~,SWIR1]=find(abs(wl-1650)==min(abs(wl-1650)));
    [~,SWIR2]=find(abs(wl-2220)==min(abs(wl-2220)));
    
    if sum(indx==1) % Clay Minerals (CM)
        iter=iter+1;
        
        IM=squeeze(M(:,:,SWIR1(1)))./squeeze(M(:,:,SWIR2(1)));
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
end

if sum(indx==2) % Sand moisture index (SMI)
        iter=iter+1;
        wlr=[1860 1925 2140];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=sqrt((squeeze(M(:,:,ind1(1)).^2+...
squeeze(M(:,:,ind2(1)).*(squeeze(M(:,:,ind3(1)).^3));
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end

if sum(indx==3) % SWIR Fine particles Index (FI)
        iter=iter+1;
        wlr=[2133 2225 2209];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind2(1)).^2/...
(squeeze(M(:,:,ind2(1)).^2+squeeze(M(:,:,ind2(1)).^2);
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
        
        for i=1:size(C,3)
        ha=[];
        IM=squeeze(C(:,:,i));
        figure;
        ha(1)=subplot(3,1,1);
        imagesc(d,d(1:size(C,1)),IM)
        caxis([nanmean(IM(:))-3*nanstd(IM(:)) nanmean(IM(:))+3*nanstd(IM(:))])
        colorbar
        map=associated_swir_map{indx(i)};
        colormap(map)
        xlabel('Depth (cm)')
        ylabel('Width (cm)')
        set(gca,'fontsize',14)
        title(index_swir{indx(i)})
        ha(2)=subplot(3,1,2);
        imagesc(d,d(1:size(IM,1)),RGB)
        colorbar
        xlabel('Depth (cm)')
        ylabel('Width (cm)')
        set(gca,'fontsize',14)
        ha(3)=subplot(313);
        plot(d,nanmean(IM(round(0.4*size(IM,1):0.6*size(IM,1)),:),1),'color',associated_swir_p{indx(i)})
        grid on
        xlabel('Depth (cm)')
        ylabel('SI')
        set(gca,'fontsize',14)
        colorbar
        xlim([min(d) max(d)])
        linkaxes(ha,'x')
    end

    
    figure;
    imagesc(corr(reshape(C,[],size(C,3))))
    colorbar
    colormap(jet)
    xtickangle(45)
    set(gca,'fontsize',14,'xtick',1:length(label),'xticklabel',label,'ytick',1:length(label),'yticklabel',label)
end
end