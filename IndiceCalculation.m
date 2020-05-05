function [C,label] = IndiceCalculation(M,RGB,d,wl,indx)
% Fonction qui permet de calculer les indices de la littérature
% INPUT :
%           M : cube hyperspectral
%           wl : longueur d'onde
% OUTPUT :
%           C : carte des indices
%           label : indices calculés

if length(wl)==98||length(wl)==144
    wl=wl(10:end-10);
    M=M(:,:,10:end-10);
end

if mean(RGB(:))>1000
    RGB=RGB/10000;
end
RGB=RGB*(0.5/mean(RGB(:)));

map_red=cat(2,[0.5:0.001:1;zeros(1,501);zeros(1,501)],[ones(1,1001);0:0.001:1;0:0.001:1])';
map_green=cat(2,[zeros(1,501);0.5:0.001:1;zeros(1,501)],[0:0.001:1;ones(1,1001);0:0.001:1])';
map_orange=[ones(1,1001);1:-0.0005:.5;zeros(1,1001)]';
map_blue=cat(2,[zeros(1,501);zeros(1,501);0.5:0.001:1],[0:0.001:1;0:0.001:1;ones(1,1001)])';
map_gray=[0:0.001:1;0:0.001:1;0:0.001:1]';
p_r=[1 0 0];
p_g=[0 1 0];
p_b=[0 0 1];
p_o=[1 1 0];
p_gr=[0.5 0.5 0.5];

M=M/10000*100; % % Reflectance

%% Chlorophyll
% Check the wavelengths range
iter=0;
if sum(wl<800&wl>600)>0
    
    index_vnir={'Chlorophylls R675/R750',...
        'Chlorophylls R645/R675',...
        'Chlorophylls R660/R670',...
        'Chlorophylls R590/R690',...
        'Chlorophylls R675-R750',...
        'Chlorophylls R645-R675',...
        'Chlorophylls dR660-dR690',...
        'Chlorophylls d675',...
        'Chlorophylls d690',...
        'Chlorophylls RABD670',...
        'Chlorophylls I-band',...
        'Chlorophylls Area650-750',...
        'Chlorophylls Area650-700',...
        'Chlorophylls Area690-730',...
        'Chlorophylls Area600-760',...
        'Chlorophylls Area600-650',...
        'Chlorophylls Area690-730',...
        'Chlorophylls (Area650-700)/R670',...
        'Phycocyanin (Area600-630)/R615',...
        'Bacteriochlorophyll a RABD846',...
        'Carotenoids RABD510',...
        'Oxydes d555',...
        'Oxydes RABA400-560',...
        'Iron Oxydes R720/R880'};
    associated_vnir_map={map_green;	map_green;	flip(map_green);	map_green;	map_green;	map_green;	map_green;	flip(map_green);	flip(map_green);	map_green;	map_green;	flip(map_green);	flip(map_green);	flip(map_green);	flip(map_green);	flip(map_green);	flip(map_green);	map_green;	map_green;	flip(map_green);	map_orange;	map_red;	map_red;	flip(map_red)};
    associated_vnir_p={p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_o;	p_r;	p_r;	p_r};
    
    if nargin<5
        [indx,~] = listdlg('ListString',index_vnir);
    end
    
    for i=1:length(indx)
        index_vnir_selec{i}=index_vnir{indx(i)};
    end
    
    if sum(double(strcmp('Chlorophylls R675/R750',index_vnir_selec)))
        %Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[675 750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls R645/R675',index_vnir_selec)))
        %Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[645 675];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls R660/R670',index_vnir_selec)))
        %von Gunten, L., Grosjean, M. (2009) High-resolution, quantitative climate reconstruction over the past 1000 years and pollution history derived from lake sediments in Central Chile. Philosophisch-naturwissenschaftlichen Fakultät 246
        iter=iter+1;
        wlr=[660 670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls R590/R690',index_vnir_selec)))
        %Trachsel, M., Grosjean, M., Schnyder, D., Kamenik, C., Rein, B. (2010) Scanning reflectance spectroscopy (380–730 nm): a novel method for quantitative high-resolution climate reconstructions from minerogenic lake sediments. Journal of Paleolimnology 44: 979–994
        iter=iter+1;
        wlr=[590 690];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls R675-R750',index_vnir_selec)))
        %Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[675 750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls R645-R675',index_vnir_selec)))
        % Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[645 675];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls dR660-dR690',index_vnir_selec)))
        % Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[660 690];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        Md=reshape(savgol(reshape(M,[],size(M,3)),7,2,1),size(M,1),size(M,2),size(M,3));
        
        IM=squeeze(Md(:,:,ind1(1)))-squeeze(Md(:,:,ind2(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls d675',index_vnir_selec)))
        % Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=675;
        Md=reshape(savgol(reshape(M,[],size(M,3)),7,2,1),size(M,1),size(M,2),size(M,3));
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        
        IM=squeeze(Md(:,:,ind1(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls d690',index_vnir_selec)))
        %Wolfe, A.P., Vinebrooke, R.D., Michelutti, N., Rivard, B., Das, B. (2006) Experimental calibration of lake-sediment spectral reflectance to chlorophyll a concentrations: methodology and paleolimnological validation. Journal of Paleolimnology 36: 91–100
        iter=iter+1;
        wlr=690;
        Md=reshape(savgol(reshape(M,[],size(M,3)),7,2,1),size(M,1),size(M,2),size(M,3));
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        
        IM=squeeze(Md(:,:,ind1(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls RABD670',index_vnir_selec)))
        % Rein, B., Sirocko, F. (2002) In-situ reflectance spectroscopy - analysing techniques for high-resolution pigment logging in sediment cores. International Journal of Earth Sciences 91: 950–954
        iter=iter+1;
        wlr=[590 730 670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(((ind3(1)-ind1(1))*squeeze(M(:,:,ind1(1)))+...
            (ind2(1)-ind3(1))*squeeze(M(:,:,ind2(1))))/...
            ((ind3(1)-ind1(1))+(ind2(1)-ind3(1))))./squeeze(M(:,:,ind3(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls I-band',index_vnir_selec)))
        % Rein, B., Sirocko, F. (2002) In-situ reflectance spectroscopy - analysing techniques for high-resolution pigment logging in sediment cores. International Journal of Earth Sciences 91: 950–954
        iter=iter+1;
        wlr=[590 730 670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(((ind3(1)-ind1(1))*squeeze(M(:,:,ind1(1)))+(ind2(1)-ind3(1))*squeeze(M(:,:,ind2(1))))/((ind3(1)-ind1(1))+(ind2(1)-ind3(1))))./squeeze(M(:,:,ind3(1)))./(squeeze(mean(M,3)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls Area650-750',index_vnir_selec)))
        % Wolfe, A.P., Vinebrooke, R.D., Michelutti, N., Rivard, B., Das, B. (2006) Experimental calibration of lake-sediment spectral reflectance to chlorophyll a concentrations: methodology and paleolimnological validation. Journal of Paleolimnology 36: 91–100
        iter=iter+1;
        wlr=[650 750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))/100*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))/100*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1))/100,3)));%*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls Area650-700',index_vnir_selec)))
        % Wolfe, A.P., Vinebrooke, R.D., Michelutti, N., Rivard, B., Das, B. (2006) Experimental calibration of lake-sediment spectral reflectance to chlorophyll a concentrations: methodology and paleolimnological validation. Journal of Paleolimnology 36: 91–100
        iter=iter+1;
        wlr=[650 700];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))/100*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))/100*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1))/100,3)));%*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls Area690-730',index_vnir_selec)))
        % Trachsel, M., Grosjean, M., Schnyder, D., Kamenik, C., Rein, B. (2010) Scanning reflectance spectroscopy (380–730 nm): a novel method for quantitative high-resolution climate reconstructions from minerogenic lake sediments. Journal of Paleolimnology 44: 979–994
        iter=iter+1;
        wlr=[590 730];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=141*squeeze(M(:,:,ind2(1)))+...
            141*(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))/2-...
            6*squeeze(M(:,:,ind1(1)))+5*squeeze(M(:,:,ind2(1)))-...
            10*sum(squeeze(M(:,:,ind1(1):ind2(1))),3)./mean(M,3);
        %         IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))/100*(wl(ind2(1))-wl(ind1(1))))/2+...
        %             squeeze(M(:,:,ind1(1)))/100*(wl(ind2(1))-wl(ind1(1)))-...
        %             squeeze(sum(M(:,:,ind1(1):ind2(1))/100,3)));%*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls Area600-760',index_vnir_selec)))
        % Das, B., Vinebrooke, R.D., Sanchez-azofeifa, A., Rivard, B., Wolfe, A.P. (2005) Inferring sedimentary chlorophyll concentrations with reflectance spectroscopy : a novel approach to reconstructing historical changes in the trophic status of mountain lakes. Canadian Journal of Fisheries and Aquatic Sciences 62: 1067–1078
        iter=iter+1;
        wlr=[600 760];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))/100*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))/100*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1))/100,3)));%*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls Area600-650',index_vnir_selec)))
        % Van Exem erreur? Modif Wolfe
        iter=iter+1;
        wlr=[600 650];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))/100*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))/100*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1))/100,3)));%*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls Area690-730',index_vnir_selec)))
        % Van Exem erreur? Modif Trashsel
        iter=iter+1;
        wlr=[690 730];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))/100*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))/100*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1))/100,3)));%*mean(wl(2:end)-wl(1:end-1));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Chlorophylls (Area650-700)/R670',index_vnir_selec)))
        % Van Exem
        iter=iter+1;
        wlr=[650 700 670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1)),3)))./(squeeze(M(:,:,ind3(1)))*100);%*mean(wl(2:end)-wl(1:end-1))
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Phycocyanin (Area600-630)/R615',index_vnir_selec)))
        iter=iter+1;
        wlr=[600 630 615];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))/100*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))/100*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1))/100,3)))./(squeeze(M(:,:,ind3(1)))/100);%*mean(wl(2:end)-wl(1:end-1))
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Bacteriochlorophyll a RABD846',index_vnir_selec)))
        iter=iter+1;
        wlr=[800 900 844]; % Butz et al., 2015
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=((55*squeeze(M(:,:,ind1(1)))+55*squeeze(M(:,:,ind2(1))))/(55+55))./squeeze(M(:,:,ind3(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Carotenoids RABD510',index_vnir_selec)))
        % Rein, B., Sirocko, F. (2002) In-situ reflectance spectroscopy - analysing techniques for high-resolution pigment logging in sediment cores. International Journal of Earth Sciences 91: 950–954
        iter=iter+1;
        wlr=[490 530 510];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(((ind2-ind3)*squeeze(M(:,:,ind1(1)))+...
            (ind3-ind1)*squeeze(M(:,:,ind2(1))))/2)./squeeze(M(:,:,ind3(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Oxydes d555',index_vnir_selec)))
        iter=iter+1;
        wlr=555;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        Md=reshape(savgol(reshape(M,[],size(M,3)),7,2,1),size(M,1),size(M,2),size(M,3));
        IM=squeeze(Md(:,:,ind1(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Oxydes RABA400-560',index_vnir_selec)))
        iter=iter+1;
        wlr=[400 560 590];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=sum(repmat(squeeze(M(:,:,ind3(1))),1,1,ind2(1)-ind1(1)+1)./M(:,:,ind1(1):ind2(1)),3)./mean(M,3);
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(double(strcmp('Iron Oxydes R720/R880',index_vnir_selec)))
        % Jackisch, R., Lorenz, S., Zimmermann, R., Möckel, R., Gloaguen, R. (2018) Drone-borne hyperspectral monitoring of acid mine drainage: An example from the Sokolov lignite district. Remote Sensing 10: 1–23
        iter=iter+1;
        wlr=[720 880];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(isnan(IM)==1);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=nanmedian(nanmedian(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
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
    
    if sum(wl<2400&wl>1100)==0
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
        imagesc(corr(reshape(C,[],size(C,3)), 'rows','complete'))
        colorbar
        colormap(jet)
        xtickangle(45)
        set(gca,'fontsize',14,'xtick',1:length(label),'xticklabel',label,'ytick',1:length(label),'yticklabel',label)
    else
        labelvnir=label;
    end
end

if sum(wl<2400&wl>1100)>0
    
    index_swir={'Clay R2200',...
        'Clay RABA2190-2230',...
        'Hydroxyl bonds R1450',...
        'Moisture R1935',...
        'R1935/R1450',...
        'Normalized Differenced Gypsum Ratio (NDGI)',...
        'Calcite R2340',...
        'Illite-muscovite Crystallinity R1900/R2200',...
        'Kaolinite Crystallinity Index (KCI)'};
    associated_swir_map={flip(map_orange);	map_orange;	map_blue;	map_blue;	map_blue;	jet(1500);	jet(1500);	jet(1500);	jet(1500)};
    associated_swir_p={p_o;	p_o;	p_b;	p_b;	p_b;	[0 0 0];	[0 0 0];	[0 0 0];	[0 0 0]};
    
    if nargin<5
        [indx,~] = listdlg('ListString',index_swir);
    end
    
    for i=1:length(indx)
        index_swir_selec{i}=index_swir{indx(i)};
    end
    
    if sum(double(strcmp('Clay R2200',index_swir_selec)))
        iter=iter+1;
        wlr=2200;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        IM=squeeze(M(:,:,ind1(1)));
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
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Clay RABA2190-2230',index_swir_selec)))
        iter=iter+1;
        wlr=[2190 2230];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(((squeeze(M(:,:,ind2(1)))-squeeze(M(:,:,ind1(1))))*(wl(ind2(1))-wl(ind1(1))))/2+...
            squeeze(M(:,:,ind1(1)))*(wl(ind2(1))-wl(ind1(1)))-...
            squeeze(sum(M(:,:,ind1(1):ind2(1)),3)))*mean(wl(2:end)-wl(1:end-1));
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
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Hydroxyl bonds R1450',index_swir_selec)))
        iter=iter+1;
        wlr=1450;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        IM=squeeze(M(:,:,ind1(1)));
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
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Moisture R1935',index_swir_selec)))
        iter=iter+1;
        wlr=1935;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        IM=squeeze(M(:,:,ind1(1)));
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
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('R1935/R1450',index_swir_selec)))
        iter=iter+1;
        wlr=[1935 1450];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
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
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Normalized Differenced Gypsum Ratio (NDGI)',index_swir_selec)))
        % Milewski, R., Chabrillat, S., Brell, M., Schleicher, A.M., Guanter, L. (2019) Assessment of the 1.75 µm absorption feature for gypsum estimation using laboratory, air- and spaceborne hyperspectral sensors. International Journal of Applied Earth Observation and Geoinformation 77: 69–83
        iter=iter+1;
        wlr=[1690 1750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
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
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Calcite R2340',index_swir_selec)))
        % Sun, L., Khan, S., Godet, A. (2018) Integrated ground-based hyperspectral imaging and geochemical study of the Eagle Ford Group in West Texas. Sedimentary Geology 363: 34–47
        iter=iter+1;
        wlr=2340;
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        IM=squeeze(M(:,:,ind1(1)));
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
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Illite-muscovite Crystallinity R1900/R2200',index_swir_selec)))
        % Doublier, M.P., Roache, T., Potel, S. (2010) Short-wavelength infrared spectroscopy: A new petrological tool in low-grade to very low-grade pelites. Geology 38: 1031–1034
        iter=iter+1;
        wlr=[2200 1900];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
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
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Kaolinite Crystallinity Index (KCI)',index_swir_selec)))
        % Alonso de Linaje, V., Khan, S.D. (2017) Mapping of diagenetic processes in sandstones using imaging spectroscopy: A case study of the Utrillas Formation, Burgos, Spain. Sedimentary Geology 353: 114–124
        iter=iter+1;
        wlr=[2120 2430 2159 2165];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [Mcr,~,wlcr]=continuum_removal(wl(ind1(1):ind2(1)),M(:,:,ind1(1):ind2(1)));
        wlcr=wlcr';
        [~,ind3]=find(abs(wlcr-wlr(1,3))==min(abs(wlcr-wlr(1,3))));
        [~,ind4]=find(abs(wlcr-wlr(1,4))==min(abs(wlcr-wlr(1,4))));
        if ind4==ind3
            ind4=ind4+1;
        end
        IM=squeeze(Mcr(:,:,ind3(1)))./squeeze(Mcr(:,:,ind4(1)));
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
    
    if sum(wl<800&wl>600)>0
        labeltmp=label;
        label=labelvnir;
        for i=1:length(labeltmp)
            label{i+length(labelvnir)}=labeltmp{i};
        end
    end
    
    figure;
    imagesc(corr(reshape(C,[],size(C,3)), 'rows','complete'))
    colorbar
    colormap(jet)
    xtickangle(45)
    set(gca,'fontsize',14,'xtick',1:length(label),'xticklabel',label,'ytick',1:length(label),'yticklabel',label)
end
end