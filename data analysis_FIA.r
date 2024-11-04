###############################################################
#### get fia data ----
###############################################################
#download FIA plot
devtools::install_github('hunter-stanke/rFIA')
state.list=as.character(read.csv("D:/Functional diversity and eco funtion/state.list.csv")[,1])
library(rFIA);require(dplyr)	
err=c()
for (i in state.list){
	#dir.create(paste0("data analysis2/FIA20240327/",i))	
	skip_to_next <- FALSE  
	tryCatch({
		getFIA(states=i, dir =paste0("data analysis2/FIA20240327/",i), common = FALSE, tables = c('PLOT', 'TREE','SURVEY','COND','TREE_GRM_COMPONENT'),load = FALSE, nCores = 7)
		#getFIA(states=i, dir =paste0("data analysis2/FIA20240327/",i), common = FALSE, tables = c('TREE_GRM_COMPONENT'),load = FALSE, nCores = 7)		
	}, error=function(e){skip_to_next <<- TRUE})
	if(skip_to_next) err=c(err,i)
}

## get abundance data
library(dplyr);require(parallel)
get.fia2=function(statename){
	require(rFIA)
	require(magrittr)
	require(dplyr)
	#state <-  readFIA(paste('data analysis2/FIA20240327/',statename,sep=""),tables=c('PLOT', 'TREE','SURVEY','COND','TREE_GRM_COMPONENT'))
	state <-  readFIA(paste('data analysis2/FIA20240327/',statename,sep=""),tables=c('PLOT', 'TREE','SURVEY','COND'))
	state$TREE$pltID <-  stringr::str_c(state$TREE$UNITCD, state$TREE$STATECD, state$TREE$COUNTYCD, state$TREE$PLOT,state$TREE$MEASYEAR,sep = '_')
	state$TREE$pltID.no.year <-  stringr::str_c(state$TREE$UNITCD, state$TREE$STATECD, state$TREE$COUNTYCD, state$TREE$PLOT,sep = '_')
	state$TREE$Abund_weight <-  (state$TREE$DIA*0.0254)^2*(state$TREE$TPA_UNADJ * 2.4710538147)   ##TPA base level TPA  in m2/ha	
	cond=select(state$COND,PLT_CN,STDORGCD,STDSZCD,COND_STATUS_CD,CONDPROP_UNADJ,PHYSCLCD,DSTRBCD1,STDAGE)
	cond$DSTRBCD1=ifelse(cond$DSTRBCD1==0,0,1)	
	cond$PHYSCLCD=ifelse(cond$PHYSCLCD<20,"Xeric",ifelse(cond$PHYSCLCD<30,"Mesic","Hydric"))
	cond=cond%>%distinct
	# state$TREE_GRM_COMPONENT=state$TREE_GRM_COMPONENT%>%select(TRE_CN,MICR_COMPONENT_AL_FOREST,SUBP_COMPONENT_AL_FOREST)
	# gc();
	tree_df <- state$TREE %>%filter(INVYR!=9999)%>%#filter(INVYR!=9999&STATUSCD==1)%>%
		left_join(select(state$PLOT,MEASYEAR,MEASMON,LAT,LON,ELEV,CN,PLOT_STATUS_CD,SRV_CN,ECOSUBCD,PREV_PLT_CN), by=c("PLT_CN"="CN")) %>% filter(PLOT_STATUS_CD ==1&(MEASYEAR!=9999)) %>%
		left_join(select(state$SURVEY,CN,ANN_INVENTORY,P3_OZONE_IND), by=c("SRV_CN"="CN")) %>% filter(ANN_INVENTORY == "Y" & P3_OZONE_IND == "N") %>% 
		#left_join(state$TREE_GRM_COMPONENT, by=c("CN"="TRE_CN")) %>% 
		left_join(cond, by="PLT_CN",relationship =  "many-to-many") %>% filter(STDORGCD == 0&STDSZCD!=5)%>%#&CONDPROP_UNADJ >= 0.95&COND_STATUS_CD==1)%>%
		select(pltID,pltID.no.year,CN,PREV_TRE_CN,SPCD,Abund_weight,MEASYEAR,MEASMON,LAT,LON,ELEV,HT,TPA_UNADJ,ECOSUBCD,STDAGE,PHYSCLCD,DRYBIO_AG,DRYBIO_BG,DSTRBCD1,PLT_CN,PREV_PLT_CN,
			RECONCILECD,DIA,CONDPROP_UNADJ,STATUSCD,COND_STATUS_CD)%>%
		filter(!is.na(Abund_weight))%>%distinct
	tree_df$ECOSUBCD=substr(tree_df$ECOSUBCD,1,nchar(tree_df$ECOSUBCD)-2)%>%gsub("^\\s|\\s$", "", .)	
	tree_df$statename=statename
	return(tree_df)#CN is an unique tree record
	rm(state,cond,tree_df);gc()
}
state.list=as.character(read.csv("D:/Functional diversity and eco funtion/state.list.csv")[,1])
# fia0=c()
# for (statename in state.list) {
	# tmp=get.fia2(statename);
	# fia0=rbind(fia0,tmp);
	# rm(tmp);gc()
	# print(which(state.list==statename))
# }
no_cores <- detectCores() - 1
mycl <- makePSOCKcluster(no_cores); 
fia0=do.call(rbind,parLapply(cl=mycl,state.list,get.fia2))
stopCluster(mycl)
sp=read.csv("data analysis2/FIA20240327/REF_SPECIES.csv")
fia=fia0%>%left_join(sp[,c("SPCD","SCIENTIFIC_NAME")],by="SPCD")
spname=unique(fia$SCIENTIFIC_NAME)

#nomenclature
library(TNRS)
res=TNRS(taxonomic_names = spname,sources = "wfo")%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1&Overall_score>=0.9)%>%
	dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
res2=read.csv("unmat.csv")%>% dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)%>%rbind(res)%>%.[order(.$Taxonomic_status),]%>% .[!duplicated(.[,'Name_submitted']),]
fia=fia%>%left_join(select(res2,Name_submitted,Accepted_name,Accepted_name_id),by=c("SCIENTIFIC_NAME"="Name_submitted"))%>%filter(!is.na(Accepted_name))

#invasive plants
res=read.csv("data analysis2/res.US-RIIS.csv")%>%select(-Name_submitted)
fia=fia%>%left_join(res,by="Accepted_name")
## convert units - conversion of per acre to per hectare
fia$TPA_UNADJ <-  fia$TPA_UNADJ*2.4710538147   ##TPA base level TPA
## convert units -- lbs dry biomass to Mg (metric tons)
fia$DRYBIO_AG <- fia$DRYBIO_AG * 0.00045359237*fia$TPA_UNADJ
fia$DRYBIO_BG <- fia$DRYBIO_BG * 0.00045359237*fia$TPA_UNADJ
fia$Biomass=fia$DRYBIO_AG+fia$DRYBIO_BG
## feet to meter;inch to cm
fia$ELEV=fia$ELEV*0.3048 #feet to meter
fia$HT=fia$HT*0.3048
fia$DIA=fia$DIA*2.54
## remove unpossible values
fia[(fia$HT<0|fia$HT>120)&(!is.na(fia$HT)),"HT"]=NA
fia[fia$Accepted_name=="Salix x pendulina nothof. x salamonii","Accepted_name"]="Salix x sepulcralis"
fia$DIA2=ifelse(fia$DIA>=12.5,"tree","sapling")
#fia$INVR = substr(fia$pltID,nchar(fia$pltID)-3,nchar(fia$pltID))
save(fia,file="data analysis2/fia.abundance.rda")

allplots=fia%>%filter(LON>=-100&MEASYEAR>=2000&MEASYEAR<=2021&STATUSCD==1&CONDPROP_UNADJ ==1)
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
inv.plots=allplots%>%filter(degreeOfEstablishment%in%status[1:2])%>%select(PLT_CN)%>%distinct
inv.plots$invaded=1
allplots=allplots%>%left_join(inv.plots,by="PLT_CN")
allplots[is.na(allplots$invaded),"invaded"]=0

# plt.pre=allplots%>%select(PLT_CN,PREV_PLT_CN)%>%distinct%>%na.omit
# allplots.pre=fia%>%filter(PLT_CN%in%plt.pre$PREV_PLT_CN)
# allplots.after=fia%>%filter(PREV_PLT_CN%in%plt.pre$PLT_CN)

# invaded.plots=fia%>%filter(LON>=-100&invaded==1)
# plt.pre=invaded.plots%>%select(PLT_CN,PREV_PLT_CN)%>%distinct%>%na.omit
# invaded.pre=fia%>%filter(PLT_CN%in%plt.pre$PREV_PLT_CN)
# invaded.after=fia%>%filter(PREV_PLT_CN%in%plt.pre$PLT_CN)

# allplots=rbind(allplots.pre,allplots,allplots.after,invaded.plots,invaded.pre,invaded.after)%>%
	# distinct%>%filter(!is.na(Abund_weight))
save(allplots,file="data analysis2/fia.abundance.ena2.rda")

##statistics
load("data analysis2/fia.abundance.rda")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=fia%>%filter(STATUSCD%in%c(1,2)&LON>=-100&CONDPROP_UNADJ==1&(!degreeOfEstablishment%in%status[3]))#STATUSCD==1
allplots%>%select(PLT_CN)%>%distinct%>%nrow
allplots%>%filter(!is.na(degreeOfEstablishment))%>%select(Accepted_name,pltID)%>%distinct%>%dim
allplots%>%filter(!is.na(degreeOfEstablishment))%>%select(Accepted_name,pltID)%>%distinct%>%group_by(Accepted_name)%>%tally%>%.[order(.$n,decreasing=T),]

#########################################
### characters of invasive speices based on USDA plants,such as whether shade tolerance, height, where are these invasive species come from and how long ----
#########################################
library(dplyr);library(sf);library(ggplot2)
ecoregion=st_read("S_USA.EcoMapProvinces/S_USA.EcoMapProvinces.shp")%>% st_transform("+proj=longlat +datum=WGS84")
load("data analysis2/fia.abundance.ena2.rda")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
inv.splist=allplots%>%filter(degreeOfEstablishment%in%status[1:2])%>%select(Accepted_name)%>%distinct
#write.csv(inv.splist,"data analysis2/inv.splist.csv")#then mannual add characters of species

splist=inv.splist$Accepted_name
theme=theme_bw()+
		theme(axis.text = element_text(size=12,color='black'),			
			axis.title = element_text(size=15,color='black'),
			axis.text.x=element_text(angle=30),
			strip.text=element_text(size=12,color='black',face="bold.italic"),
		strip.background = element_rect(fill = 'white', colour = 'black', linewidth = rel(2), linetype = 2))	

#plot number of occupied plots changes through time
all.inv=allplots%>%filter(Accepted_name%in%splist)%>%group_by(MEASYEAR,Accepted_name)%>%tally
ggplot(all.inv, aes(x=MEASYEAR, y=n)) + 	
		geom_point(size=3,shape=21,color="black",fill="#F8766D",show.legend=FALSE) + 		
		labs(x="Year",y="Plots occupied") +
		#geom_smooth(aes(x = MEASYEAR, y = n),data=all.inv,method="gam",linewidth=1,show.legend=FALSE,se =F,linetype=1)+		
		facet_wrap(~Accepted_name, scales="fixed",nrow=4)+
		theme
		
#map each invasive species
plots.sp=allplots%>%filter(Accepted_name%in%splist)%>% dplyr::select(pltID.no.year,LON,LAT,Accepted_name)%>% distinct%>%
	st_as_sf(coords=c("LON","LAT"))%>% st_set_crs("+proj=longlat +datum=WGS84")
ggplot() + labs(x="Longitude",y="Latitude")+		
		geom_sf(data=ecoregion,color="lightgray",fill ="white") +xlim(-100,-67)+ 
		geom_sf(data=plots.sp,color="#F8766D",size=1,show.legend=F)+
		facet_wrap(~Accepted_name, scales="fixed",nrow=4)+
		theme+theme(panel.background = element_rect(fill = '#619CFF', colour = 'black'),axis.text=element_text(size=10))

#########################################
##### map plots and invaded plots of ENA  #########
#########################################
library(dplyr);library(ggpubr);library(sf);
load("data analysis2/fia.abundance.ena.rda")
plots.sp=allplots%>%filter(LON>=-100&MEASYEAR>=2000&MEASYEAR<=2021) %>% dplyr::select(pltID.no.year,LON,LAT,invaded)%>% distinct%>%
	st_as_sf(coords=c("LON","LAT"))%>% st_set_crs("+proj=longlat +datum=WGS84")
plots.sp$invaded2=ifelse(plots.sp$invaded==0,"Uninvaded","Invaded")
plots.sp$invaded2=factor(plots.sp$invaded2,levels=c("Uninvaded","Invaded"))

#part1 study area
US=map_data("world", region = "USA")%>%filter(!subregion%in%c("Hawaii","Alaska"))
eUS=US%>%filter(long>-100)
p=ggplot() +  
	geom_polygon(data=US,aes(x=long,y=lat,group=group),color = NA,fill="lightgray",size = 1) + 
	geom_polygon(data=eUS,aes(x=long,y=lat,group=group),color = NA,fill="black",size = 1) +
	coord_equal()+theme_bw()+ theme(axis.text = element_blank(),axis.title = element_blank(),
	panel.grid.major =element_blank(),axis.ticks=element_blank(),
	panel.grid.minor =element_blank())+  geom_path()+ coord_map("albers",lat0=39, lat1=45)
	
#part 2 main map
ecomap=st_read("S_USA.EcoMapProvinces/S_USA.EcoMapProvinces.shp")%>% st_transform("+proj=longlat +datum=WGS84")
the=theme_bw()+theme(panel.background = element_rect(fill = NA, colour = 'black'),
			axis.text = element_text(size=12,color='black'),
			axis.title = element_blank(),
			panel.grid =element_blank(),
			legend.position = c(0.85,0.3),
			legend.background=element_rect(fill='transparent'),
			legend.text=element_text(face="bold",size=12),
			legend.title=element_blank(),
			legend.key = element_blank())
ggplot() +  
		geom_sf(data=plots.sp,size=0.8,alpha=1,aes(color=invaded2))+
		geom_sf(data=ecomap,color="black",fill="transparent",show.legend=F) +				
		xlim(-100,-67)+
		labs(color="Type of plot")+
		scale_color_manual(values=c("#619CFF","#F8766D"))+
		guides(color=guide_legend(ncol=1,byrow=T,override.aes = list(size = 3)))+	
		the+patchwork::inset_element(p, left = 0.61, bottom = 0.01, right = 0.95, top = 0.2)

## plots for each year	----
library(dplyr)
load("data analysis2/fia.abundance.ena.rda")
#ecoregion ref:Ecological Subregions: Sections and Subsections for the Conterminous United States https://data.fs.usda.gov/geodata/edw/datasets.php?dsetParent=EcomapSections_2007
ecoregion=st_read("S_USA.EcoMapProvinces/S_USA.EcoMapProvinces.shp")%>% st_transform("+proj=longlat +datum=WGS84")
gwr.map=c()
for(yr in unique(allplots$MEASYEAR)){
	tmp <- allplots%>%filter(MEASYEAR==yr) %>% dplyr::select(pltID.no.year,LON,LAT,invaded,MEASYEAR)%>% distinct
	gwr.map=rbind(gwr.map,tmp)
}	
plt.inv=gwr.map%>%st_as_sf(coords=c("LON","LAT"))%>% st_set_crs("+proj=longlat +datum=WGS84")
plt.inv$invaded=factor(plt.inv$invaded,levels=c(1,0))	
ggplot() + geom_sf(data=plt.inv%>%filter(invaded==0),aes(color=invaded),size=1,alpha=1)+
		geom_sf(data=plt.inv%>%filter(invaded==1),aes(color=invaded),size=1,alpha=1)+
		labs(x="Longitude",y="Latitude",color="Type of plot")+
		geom_sf(data=ecoregion,color="black",fill = NA,alpha=0.8) +xlim(-100,-67)+ 
		scale_color_manual(values=c("0"="#619CFF", "1"="#F8766D"),labels=c("Uninvaded","Invaded"))+
		#annotate("text",x=-72 , y= 27,size=2.5,	label = yr,color = "black")+
		facet_wrap(~MEASYEAR)+
		theme_bw()+
		guides(color = guide_legend(override.aes = list(size = 3))) +
		theme(axis.text = element_text(size=8,color='black'),			
			axis.title = element_text(size=15,color='black'),
			axis.text.x=element_text(angle=30),
			strip.text=element_text(size=12,color='black',face="bold.italic"),
		strip.background = element_rect(fill = 'white', colour = 'black', linewidth = rel(2), linetype = 2),
		legend.position=c(0.91,0.1),
		legend.background = element_rect(fill = NA),
		legend.text=element_text(face="bold",size=12),
		legend.title=element_blank())	
		
######################################################
##richness and biomass changes through time #######
######################################################
library(dplyr);library(ggpubr);require(ggpmisc)
#detach(package:raster)
load("data analysis2/fia.abundance.ena2.rda")
# allplots=allplots%>%filter(LON>=-100&MEASYEAR>=2000&MEASYEAR<=2021)

get.stat=function(allplots.t,vars,ecoregion=FALSE,status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")){
	colnames(allplots.t)[colnames(allplots.t)==vars]="VARS"	
	allplots.t=allplots.t%>%filter(!is.na(VARS))
	if (ecoregion){
		if (vars=="Accepted_name"){
			inv.plots.abd=allplots.t%>%filter(invaded==1)%>%group_by(MEASYEAR,pltID,ECOSUBCD)%>%
					summarize(nat.abd=length(unique(VARS[which(is.na(degreeOfEstablishment))])),inv.abd=length(unique(VARS[which(degreeOfEstablishment%in%status[1:2])])),.groups="keep")
			uninv.plots.abd=allplots.t%>%filter(invaded==0&is.na(degreeOfEstablishment))%>%distinct%>%group_by(MEASYEAR,pltID,ECOSUBCD)%>%summarize(n=length(unique(VARS)),.groups="keep")
		}else{
			inv.plots.abd=allplots.t%>%filter(invaded==1)%>%group_by(MEASYEAR,pltID,ECOSUBCD)%>%
					summarize(nat.abd=sum(VARS[which(is.na(degreeOfEstablishment))]),inv.abd=sum(VARS[which(degreeOfEstablishment%in%status[1:2])]),.groups="keep")
			uninv.plots.abd=allplots.t%>%filter(invaded==0&is.na(degreeOfEstablishment))%>%distinct%>%group_by(MEASYEAR,pltID,ECOSUBCD)%>%summarize(n=sum(VARS),.groups="keep")	
		}
		abd.inv=inv.plots.abd%>%group_by(MEASYEAR,ECOSUBCD)%>%
			summarize(nPlots=length(pltID),native.mean=mean(nat.abd),inv.mean=mean(inv.abd),native.se=sd(nat.abd)/sqrt(length(nat.abd)),inv.se=sd(inv.abd)/sqrt(length(inv.abd)),.groups="keep")
		abd.nat=uninv.plots.abd%>%group_by(MEASYEAR,ECOSUBCD)%>%summarize(nPlots=length(pltID),univ.native.mean=mean(n),univ.native.se=sd(n)/sqrt(length(n)),.groups="keep")
		native=abd.inv%>%select(MEASYEAR,ECOSUBCD,nPlots,native.mean,native.se);inv=abd.inv%>%select(MEASYEAR,ECOSUBCD,nPlots,inv.mean,inv.se)
		colnames(abd.nat)=colnames(native)=colnames(inv)=c("MEASYEAR","ECOSUBCD","nPlots","abd.mean","abd.se")	
	}else{
		if (vars=="Accepted_name"){
			inv.plots.abd=allplots.t%>%filter(invaded==1)%>%group_by(MEASYEAR,PLT_CN)%>%
					summarize(nat.abd=length(unique(VARS[which(is.na(degreeOfEstablishment))])),inv.abd=length(unique(VARS[which(degreeOfEstablishment%in%status[1:2])])),.groups="keep")
			uninv.plots.abd=allplots.t%>%filter(invaded==0&is.na(degreeOfEstablishment))%>%distinct%>%group_by(MEASYEAR,PLT_CN)%>%summarize(n=length(unique(VARS)),.groups="keep")
		}else{
			inv.plots.abd=allplots.t%>%filter(invaded==1)%>%group_by(MEASYEAR,PLT_CN)%>%
					summarize(nat.abd=sum(VARS[which(is.na(degreeOfEstablishment))]),inv.abd=sum(VARS[which(degreeOfEstablishment%in%status[1:2])]),.groups="keep")
			uninv.plots.abd=allplots.t%>%filter(invaded==0&is.na(degreeOfEstablishment))%>%distinct%>%group_by(MEASYEAR,PLT_CN)%>%summarize(n=sum(VARS),.groups="keep")	
		}
		abd.inv=inv.plots.abd%>%group_by(MEASYEAR)%>%summarize(nPlots=length(PLT_CN),native.mean=mean(nat.abd),inv.mean=mean(inv.abd),native.se=sd(nat.abd)/sqrt(length(nat.abd)),inv.se=sd(inv.abd)/sqrt(length(inv.abd)))
		abd.nat=uninv.plots.abd%>%group_by(MEASYEAR)%>%summarize(nPlots=length(PLT_CN),univ.native.mean=mean(n),univ.native.se=sd(n)/sqrt(length(n)))
		native=abd.inv%>%select(MEASYEAR,nPlots,native.mean,native.se);inv=abd.inv%>%select(MEASYEAR,nPlots,inv.mean,inv.se)
		colnames(abd.nat)=colnames(native)=colnames(inv)=c("MEASYEAR","nPlots","abd.mean","abd.se")	
	}
	#mean richness of plots through time
	abd.stat=rbind(data.frame(Type0="Native",Type="Native (Invaded plots)",native),data.frame(Type0="Non-native",Type="Non-native (Invaded plots)",inv),data.frame(Type0="Native",Type="Native (Uninvaded plots)",abd.nat))
	abd.stat[is.na(abd.stat)]=0
	abd.stat$Type0=factor(abd.stat$Type0,levels=c("Non-native","Native"))
	abd.stat$Type=factor(abd.stat$Type,levels=c("Non-native (Invaded plots)","Native (Invaded plots)","Native (Uninvaded plots)"))
	return(abd.stat)
}
#ba.stat=get.stat(allplots,"Abund_weight")# basal area through time
bio.stat=get.stat(allplots,"Biomass")
sp.stat=get.stat(allplots,"Accepted_name")
dens.stat=get.stat(allplots,"TPA_UNADJ")
theme=theme_bw()+
	theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),		
		strip.text = element_text(colour = 'black', size = rel(1.5)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))
get.p=function(d){
	pvals = do.call(rbind,lapply(levels(d$Type0), function(i) {
	  do.call(rbind,lapply(levels(d$Type), function(j) {   
		  if(nrow(d[d$Type0==i & d$Type==j, ]) > 1) {
			data.frame(Type0=i, Type=j, 
					   p.val=coef(summary(lm(abd.mean ~ MEASYEAR, data = d[d$Type0==i &d$Type==j, ])))[2,4])
		  }    
	  }))
	}))
	# Keep only rows with p.val < 0.05
	pvals = pvals[pvals$p.val < 0.05, ]
	d.subset = d%>%filter(Type0%in%pvals$Type0 & Type%in% pvals$Type) 
	return(d.subset)
}		


p.dens=ggplot(dens.stat, aes(x=MEASYEAR, y=abd.mean,fill=Type,shape=Type)) + 	
		geom_errorbar(aes(ymin=abd.mean-1.96*abd.se, ymax=abd.mean+1.96*abd.se,color=Type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
		scale_fill_manual(values=c("#F8766D","#00BA38", "#619CFF"))+	
 		scale_shape_manual(values=c(21,21,24))+ 	
		geom_point(size=3,color="black",show.legend=T) + 		
		labs(x="Year",y="Mean tree density within plot (trees/Ha)") +
		geom_smooth(aes(x = MEASYEAR, y = abd.mean,color=Type),data=get.p(dens.stat),method="lm",linewidth=1.5,show.legend=FALSE,se =T,linetype=1)+		
		#stat_poly_eq(aes(label =  paste(stat(p.value.label), "*\"\"",sep = ""),color=Type),rr.digits = 2,p.digits = 2,label.x=0.1,label.y=c(1,0.92),size=6)+
		facet_wrap(~Type0, scales="free_y",ncol = 1)+
		theme+theme(legend.position= c(0.3,0.95),
			legend.background = element_rect(fill = NA),
			legend.title=element_blank(),
			legend.text=element_text(face="bold",size=15))	
p.sp=ggplot(sp.stat, aes(x=MEASYEAR, y=abd.mean,fill=Type,shape=Type)) + 	
		geom_errorbar(aes(ymin=abd.mean-1.96*abd.se, ymax=abd.mean+1.96*abd.se,color=Type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
		scale_fill_manual(values=c("#F8766D","#00BA38", "#619CFF"))+	
 		scale_shape_manual(values=c(21,21,24))+ 	
		geom_point(size=3,color="black",show.legend=F) + 		
		labs(x="Year",y="Mean species richness within plot") +
		geom_smooth(aes(x = MEASYEAR, y = abd.mean,color=Type),data=get.p(sp.stat),method="lm",linewidth=1.5,show.legend=FALSE,se =T,linetype=1)+		
		#stat_poly_eq(aes(label =  paste(stat(p.value.label), "*\"\"",sep = ""),color=Type),rr.digits = 2,p.digits = 2,label.x=0.1,label.y=c(1,0.92),size=6)+
		facet_wrap(~Type0, scales="free_y",ncol = 1)+
		theme
p.bio=ggplot(bio.stat, aes(x=MEASYEAR, y=abd.mean*1.5,fill=Type,shape=Type)) + 	
		geom_errorbar(aes(ymin=abd.mean*1.5-1.96*abd.se, ymax=abd.mean*1.5+1.96*abd.se,color=Type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
		scale_fill_manual(values=c("#F8766D","#00BA38", "#619CFF"))+	
 		scale_shape_manual(values=c(21,21,24))+ 
		geom_point(size=3,color="black",show.legend=F) + 		
		labs(x="Year",y="Mean biomass within plot (Mg/Ha)") +
		geom_smooth(aes(x = MEASYEAR, y = abd.mean*1.5,color=Type),data=get.p(bio.stat),method="lm",linewidth=1.5,show.legend=FALSE,se =T,linetype=1)+
		#stat_poly_eq(aes(label =  paste(stat(p.value.label), "*\"\"",sep = ""),color=Type),rr.digits = 2,p.digits = 2,label.x=0.1,label.y=c(1,0.92),size=6)+			
		facet_wrap(~Type0, scales="free_y",ncol = 1)+
		theme
ggarrange(p.dens,p.sp,p.bio,nrow=1,labels="auto",common.legend=T)	
	
mf=function(div,dat.plt){
	dat=dat.plt%>%filter(Type==div)
	m=lm(abd.mean~MEASYEAR,dat)
	re=c(Type=div,df= summary(m)$df[2],F = round(summary(m)$f[1],2),p=round(summary(m)$coef[2,4],2))
	return(re)
}
do.call(rbind,lapply(c("Native (Invaded plots)","Non-native (Invaded plots)","Native (Uninvaded plots)"),mf,dens.stat))
do.call(rbind,lapply(c("Native (Invaded plots)","Non-native (Invaded plots)","Native (Uninvaded plots)"),mf,sp.stat))
do.call(rbind,lapply(c("Native (Invaded plots)","Non-native (Invaded plots)","Native (Uninvaded plots)"),mf,bio.stat))

model2 <- lm(abd.mean ~ MEASYEAR + Type + Type:MEASYEAR, data = bio.stat[bio.stat$Type0=="Native",])
summary(model2)$df[2];round(summary(model2)$f[1],2)

# density of non-natives (individuals per hectare), biomass of non-natives (Mg/ha), and each of these two variables expressed as a percent of plot-level density or biomass
stat=allplots%>%group_by(MEASYEAR,PLT_CN)%>%summarise(dens=sum(TPA_UNADJ[which(!is.na(degreeOfEstablishment))]),dens.por=dens/sum(TPA_UNADJ)*100,
	sr=length(unique(Accepted_name[which(!is.na(degreeOfEstablishment))])),sr.por=sr/length(unique(Accepted_name))*100,
	bio=sum(Biomass[which(!is.na(degreeOfEstablishment))]),bio.por=bio/sum(Biomass)*100)%>%na.omit
stat.mean=stat%>%group_by(MEASYEAR)%>%summarise(dens=mean(dens),dens.por=mean(dens.por),sr=mean(sr),sr.por=mean(sr.por),bio=mean(bio),bio.por=mean(bio.por))%>%
	tidyr::pivot_longer(cols=dens:bio.por,names_to = "index", values_to = "value")
stat.se=stat%>%group_by(MEASYEAR)%>%
	summarise(dens=sd(dens)/sqrt(length(PLT_CN)),dens.por=sd(dens.por)/sqrt(length(PLT_CN)),sr=sd(sr)/sqrt(length(PLT_CN)),sr.por=sd(sr.por)/sqrt(length(PLT_CN)),bio=sd(bio)/sqrt(length(PLT_CN)),bio.por=sd(bio.por)/sqrt(length(PLT_CN)))%>%
	tidyr::pivot_longer(cols=dens:bio.por,names_to = "index", values_to = "se")
stat2=stat.mean%>%left_join(stat.se,by=c("MEASYEAR","index"))	
stat2$index=factor(stat2$index,levels=unique(stat2$index))
lab=c(dens="Non-native Density (trees/ha)",dens.por="Percent of plot-level density (%)",bio="Non-native Biomass (Mg/ha)",bio.por="Percent of plot-level biomass (%)",
sr="Non-native richness",sr.por="Percent of plot-level richness (%)")
ggplot(stat2, aes(x=MEASYEAR, y=value)) + 	
	geom_errorbar(aes(ymin=value-1.96*se, ymax=value+1.96*se),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
	geom_point(size=3,shape=21,fill="#00BA38",color="black",show.legend=T) + 		
	labs(x="Year",y="Mean values across plots") +
	geom_smooth(aes(x = MEASYEAR, y = value),data=stat2,method="lm",linewidth=1.5,show.legend=FALSE,se =T,linetype=1)+		
	facet_wrap(~index, scales="free_y",ncol = 2,labeller = labeller(index=lab))+
	theme

mf=function(div,stat2){
	dat=stat2%>%filter(index==div)
	m=lm(value~MEASYEAR,dat)
	re=c(Type=as.character(div),df= summary(m)$df[2],F = round(summary(m)$f[1],2),p=round(summary(m)$coef[2,4],2))
	return(re)
}	
do.call(rbind,lapply(unique(stat2$index),mf,stat2))

## time trends in number of invaded plots
statn=allplots%>%select(PLT_CN,MEASYEAR)%>%distinct%>%group_by(MEASYEAR)%>%tally
stat.inv=allplots%>%filter(invaded==1)%>%select(PLT_CN,MEASYEAR)%>%distinct%>%group_by(MEASYEAR)%>%tally

stat.nplots=rbind(data.frame(type="All plots",statn),data.frame(type="Invaded plots",stat.inv))
ggplot(stat.nplots, aes(x=MEASYEAR, y=n)) +	
 	geom_point(size=3,color="black",fill="#00BA38",shape=21,show.legend=F) + 		
	labs(x="Year",y="Number of plots") +
	facet_wrap(~type, scales="free_y",ncol = 1)+
	theme

## plot Ecoregion map
library(sf)
ecomap=st_read("S_USA.EcoMapProvinces/S_USA.EcoMapProvinces.shp")%>% st_transform("+proj=longlat +datum=WGS84")%>%filter(MAP_UNIT_S!="Water")
nc3_points <- sf::st_point_on_surface(ecomap)
nc3_coords <- as.data.frame(sf::st_coordinates(nc3_points))
nc3_coords$MAP_UNIT_S <- ecomap$MAP_UNIT_S
nc3_coords[nc3_coords$MAP_UNIT_S=="331","Y"]=37
nc3_coords[nc3_coords$MAP_UNIT_S=="331","X"]=-102

#nc3_coords=nc3_coords%>%filter(MAP_UNIT_S%in%ecoregion)
p1=ggplot() + geom_sf(data=ecomap,color="black",aes(fill=MAP_UNIT_S),show.legend=F) +	
		#scale_fill_manual(values=c("white","#619CFF"))+
		geom_text(data = nc3_coords, aes(X, Y, label = MAP_UNIT_S), colour = "white",fontface = "bold",size=4.5)+
		xlim(-102,-67)+theme_bw()+
		theme(panel.background = element_rect(fill = "lightgray", colour = 'black'),
			axis.text = element_text(size=12,color='black'),
			axis.title = element_blank(),
			panel.grid =element_blank())	
#p1+patchwork::inset_element(p2, left = 0.61, bottom = 0.01, right = 0.84, top = 0.3)
#part1 study area
US=map_data("world", region = "USA")%>%filter(!subregion%in%c("Hawaii","Alaska"))
eUS=US%>%filter(long>-100)
p2=ggplot() +  
	geom_polygon(data=US,aes(x=long,y=lat,group=group),color = NA,fill="lightgray",size = 1) + 
	geom_polygon(data=eUS,aes(x=long,y=lat,group=group),color = NA,fill="black",size = 1) +
	coord_equal()+theme_bw()+ theme(axis.text = element_blank(),axis.title = element_blank(),
	panel.grid.major =element_blank(),axis.ticks=element_blank(),
	panel.grid.minor =element_blank())+  geom_path()+ coord_map("albers",lat0=39, lat1=45)
	
#part 2 main map
p1+patchwork::inset_element(p2, left = 0.7, bottom = 0.01, right = 0.99, top = 0.2)
	
##############################################################
##impute missing trait and phylogenty values #################
##############################################################
library(dplyr);
load("D:/Functional diversity and eco funtion/trait_mat.Rdata")
trait_mat=trait_mat%>%select(LMA,Le.N,WD,Hmax,RootDep.)%>%as.data.frame
trait_mat$sp=rownames(trait_mat)
load("data analysis/traits.all.RData")
trait_mat0=traitdata %>% tidyr::pivot_wider(names_from = taxon_name, values_from = StdValue) %>% tibble::column_to_rownames("trait") %>%t()%>%as.data.frame
trait_mat0$sp=rownames(trait_mat0)
colnames(trait_mat0)=colnames(trait_mat)
trait_mat=rbind(trait_mat,trait_mat0)%>%distinct
trait_mat[trait_mat$sp=="Salix x pendulina nothof. tristis","sp"]="Salix x sepulcralis"

load("data analysis2/fia.abundance.rda")
allplots=fia%>%filter(LON>=-100)
sp=allplots%>%select(SPCD,SCIENTIFIC_NAME,Accepted_name)%>%distinct%>%left_join(trait_mat,by=c("Accepted_name"="sp"),relationship ="many-to-many")
wd=read.csv("data analysis2/FIA20240327/REF_SPECIES.csv")%>% select(SPCD,WOOD_SPGR_GREENVOL_DRYWT)
colnames(wd)[2]="WD"
height=fia %>% group_by(Accepted_name) %>% summarize(Hmax = quantile(HT,0.99,na.rm=TRUE),.groups ="keep")
sp=sp%>%select(-Hmax,-WD)%>%left_join(height,by="Accepted_name")%>%left_join(wd,by="SPCD")
sp=sp%>%group_by(Accepted_name)%>%summarise(across(LMA:WD, ~ mean(.x, na.rm = TRUE)))
colnames(sp)[4]="RootDep"
save(sp,file="data analysis2/traits.rda")#4 out of spp missing in LMA and Le.N; 71 spp missing in rootdep

#phylogeny
library(ape);library(TNRS)
tre0=read.tree("data analysis2/NorthAmerica.ExaML.treePL.best.tre")#Park DS,2020. Replication code and data, Zenodo. doi: 10.5281/zenodo.3755913.
tip=data.frame(tip=sort(tre0$tip))
tip$Taxon_name=gsub("_", " ",tip$tip)	
res=TNRS(taxonomic_names = tip$Taxon_name,sources = "wfo")%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1&Overall_score>=0.9)%>%
	dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
tip=tip%>%left_join(res,by=c("Taxon_name"="Name_submitted"))	
tip[is.na(tip$Accepted_name),"Accepted_name"]=tip[is.na(tip$Accepted_name),"Taxon_name"]
#tip%>%filter(Accepted_name=="Salix x pendulina nothof. tristis")

# imputate missing species in phylogeny
library(phytools)
mf=function(i){re=ifelse (i[1]=="×",i[2],i[1]); return(re)}
tre.genus=data.frame(tip=tre0$tip.label,genus=do.call(c,lapply(strsplit(tre0$tip.label,split="_"),mf)))
sp$genus=do.call(rbind,lapply(strsplit(sp$Accepted_name,split=" "),mf))
sp.add=sp%>%filter(!Accepted_name%in%tip$Accepted_name)%>%select(Accepted_name,genus)%>%filter(genus%in%tre.genus$genus)#35 out of 40 missing spp.
tre=tre0%>%force.ultrametric()
for (i in 1:nrow(sp.add)) {tre=add.species.to.genus(tre,sp.add$Accepted_name[i], where="root");print(i)}
tip.keep=tip%>%select(tip,Taxon_name,Accepted_name)%>%
	rbind(data.frame(tip=gsub(" ", "_",sp.add$Accepted_name),Taxon_name=sp.add$Accepted_name,Accepted_name=sp.add$Accepted_name))%>%filter(Accepted_name%in%sp$Accepted_name)
tree=keep.tip(tre,tip=tip.keep$tip)
tree$tip=tree$tip.label=tip.keep[match(tree$tip.label,tip.keep$tip),"Accepted_name"]
save(tree,file="data analysis2/phylo.imputate.Rdata")

# imputate missing trait value
library(Rphylopars)
load("data analysis2/traits.rda")
sp.add=sp%>%select(Accepted_name, LMA, Le.N, RootDep)%>%filter(Accepted_name%in%tree$tip)%>%as.data.frame
colnames(sp.add)[1]="species"
trait.imp=phylopars(sp.add ,tree, pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE)
trait_mat =	trait.imp$anc_recon%>%as.data.frame
trait_mat$Accepted_name=rownames(trait_mat)
trait_mat.impute=sp%>%select(Accepted_name, WD,Hmax)%>%left_join(trait_mat,by="Accepted_name")%>%na.omit%>%
	rbind(sp%>%select(Accepted_name, WD,Hmax,LMA,Le.N,RootDep)%>%filter(!Accepted_name%in%tree$tip))%>%as.data.frame
trait_mat.impute[trait_mat.impute$RootDep<0&(!is.na(trait_mat.impute$RootDep)),"RootDep"]=NA
sp=as.data.frame(sp)
sp[is.na(sp$LMA),"LMA"]=trait_mat.impute[trait_mat.impute$Accepted_name%in%sp[is.na(sp$LMA),"Accepted_name"],"LMA"]
sp[is.na(sp$Le.N),"Le.N"]=trait_mat.impute[trait_mat.impute$Accepted_name%in%sp[is.na(sp$Le.N),"Accepted_name"],"Le.N"]
sp[is.na(sp$RootDep),"RootDep"]=trait_mat.impute[trait_mat.impute$Accepted_name%in%sp[is.na(sp$RootDep),"Accepted_name"],"RootDep"]
trait_mat=sp
rownames(trait_mat)=trait_mat$Accepted_name
trait_mat=as.matrix(trait_mat[,-c(1,7)])
save(trait_mat,file="data analysis2/traits.impute.rda")

##################################################################################
##PD, TD, MPD, MTD and eveness, FRic, etc ########################################
##################################################################################
library(data.table);library(dplyr);library(vegan);library(hillR);library(chemodiv);library(ape);library(FD)
load("data analysis2/fia.abundance.ena2.rda")
load("data analysis2/traits.impute.rda")
load("data analysis2/phylo.imputate.Rdata")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
#remove established species
fia=allplots%>%filter(LON>=-100&(is.na(degreeOfEstablishment)|degreeOfEstablishment%in%status[1:2]))
fia$invasive=ifelse(is.na(fia$degreeOfEstablishment),"native","invasive")

#make abundance matrix
abund_mat <- fia %>%  group_by(PLT_CN,Accepted_name)  %>% summarize(abund_w = sum(Abund_weight),.groups ="keep")%>% 
		tidyr::pivot_wider(names_from = Accepted_name, values_from = abund_w) %>% tibble::column_to_rownames("PLT_CN")	
abund_mat[is.na(abund_mat)]=0
comm <- as.matrix(abund_mat)    
comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/") #abd(i)/sum(abd)
inv.sp=fia%>%filter(invasive=="invasive")%>%select(Accepted_name)%>%distinct
comm.inv=comm[,inv.sp$Accepted_name];comm.inv=comm.inv[rowSums(comm.inv)>0,]
comm.nat=comm[,!colnames(comm)%in%inv.sp$Accepted_name];comm.nat=comm.nat[rowSums(comm.nat)>0,]
## (Phylogenetic) Shannon–Wiener and simpson and Hill Evenness (defined by equation 8 in Tuomisto 2012,Oikos 121: 1203-1218) ----
get.div1=function(comm,tree,trait_mat){
	shannon=diversity(comm);
	simpson=diversity(comm,index="simpson")
	phy.shannon=hill_phylo(comm, tree, q = 1);phy.shannon[is.infinite(phy.shannon)] <- NA
	phy.simpson=hill_phylo(comm, tree, q = 2)
	PD=hill_phylo(comm, tree, q = 0)#equal to picante::pd(), do not account for abundance	
	hilleven=chemodiv::calcDiv(sampleData = comm, type = "HillEven")
	rownames(hilleven)=rownames(comm)
	FRic=fundiversity::fd_fric(trait_mat,comm) # hull volume of functional space; Some sites had less species than traits so returned FRic is 'NA'
	colnames(FRic)[1]="PLT_CN"	
	div1=cbind(PD,shannon,simpson,phy.shannon,phy.simpson,hilleven,FRic)
	return(div1)
}
div1=get.div1(comm,tree,trait_mat)	
div1.inv=get.div1(comm.inv,tree,trait_mat)	
div1.nat=get.div1(comm.nat,tree,trait_mat)	
## PD, MPD, and phylogenetic Shannon–Wiener and simpson ----
#functional attribute diversity (FAD);ref:Plant Attribute Diversity, Resilience,and Ecosystem Function: The Nature and Significance of Dominant and Minor Species		
div.cal=function(k,tree,trait_mat){
	require(dplyr);require(ape);require(FD)
	df2 <- k[k > 0];
	SR=length(df2)	
	df2.phy=df2[names(df2)%in%tree$tip.label]
	if (length(df2.phy)<=1){
		M_p = M_Ap = NA
	}else{
		tre.plt=keep.tip(tree, tip=names(df2.phy))
		dis2=cophenetic(tre.plt)#pairwise branch length
		#Standardized mean pairwise distance, range (0,1)
		M_p = sum(dis2)/(nrow(dis2)*(nrow(dis2)-1)) 
		#Standardized mean pairwise distance weighted by abundance, Q/S(S-1), Q, Rao’s quadratic entropy; S, Species richness
		M_Ap = sum(df2.phy %*% dis2 %*% matrix(df2.phy, ncol = 1))/(nrow(dis2)*(nrow(dis2)-1))	
	}
	if (length(df2)<=1){
		M_t = M_At = FAD =NA
	}else{	
		trait_mat.t=trait_mat[names(df2),]
		gowdis_corect=function(n) {
			n[which.max(n)]=ifelse(max(n,na.rm = T)==min(n,na.rm = T)&length(na.omit(n))>1,
			max(n,na.rm = T)+min(n,na.rm = T)/10000,max(n,na.rm = T));
			return(n)
		}
		trait_mat.t=apply(trait_mat.t,2,gowdis_corect)
		G.trait.mat=gowdis(trait_mat.t,ord = "podani")
		if(sum(G.trait.mat)>0){
			dis2 <- as.matrix(G.trait.mat)
			#Standardized mean pairwise distance, range (0,1)
			M_t = sum(dis2)/(nrow(dis2)*(nrow(dis2)-1)) 
			#Standardized mean pairwise distance weighted by abundance, Q/S(S-1), Q, Rao’s quadratic entropy; S, Species richness
			M_At = sum(df2 %*% dis2 %*% matrix(df2, ncol = 1))/(nrow(dis2)*(nrow(dis2)-1))						
		}else{
			M_t = M_At = NA
		}
		FAD=sum(G.trait.mat)	
	}
	re=data.frame(SR, M_p, M_Ap, M_t, M_At, FAD)
	return(re)	
}
# fin=c()
# for (i in 116465:nrow(comm)){
	# k=comm[i,]
	# tmp=div.cal(k,tree,trait_mat)
	# fin=rbind(fin,tmp)
# }
library(parallel)
no_cores <-detectCores() - 1
mycl <- makePSOCKcluster(no_cores);		
fin=do.call(rbind,parApply(cl=mycl,comm,1,div.cal,tree,trait_mat))
fin.inv=do.call(rbind,parApply(cl=mycl,comm.inv,1,div.cal,tree,trait_mat))
fin.nat=do.call(rbind,parApply(cl=mycl,comm.nat,1,div.cal,tree,trait_mat))
stopCluster(mycl)
		

fin.inv=cbind(div1.inv,fin.inv);fin.nat=cbind(div1.nat,fin.nat)
fin=cbind(div1,fin)
colnames(fin.inv)[-7]=paste(colnames(fin.inv)[-7],"invasive",sep=".")
colnames(fin.nat)[-7]=paste(colnames(fin.nat)[-7],"native",sep=".")
div.all=list(fin,fin.inv,fin.nat)%>% purrr::reduce(left_join,by='PLT_CN')

mf=function(x){
	stat=length(na.omit(x))
	if (stat>0){
		if (stat==1) return(x) 
		if (stat>1) {a=table(x);return(names(a)[which(a==max(a))][1])}
	}else {return(NA)}
}
fia$PLT_CN=as.character(fia$PLT_CN)
plots.div=fia%>%group_by(PLT_CN)%>%summarize(pltID.no.year=unique(pltID.no.year),year=unique(MEASYEAR),invaded=unique(invaded),
	ecoregion=mf(ECOSUBCD),soiltype=mf(PHYSCLCD),disturb=mean(DSTRBCD1,na.rm=T),
	biomass=mean(DRYBIO_AG+DRYBIO_BG,na.rm=T),std=mean(STDAGE,na.rm=T),BA=sum(Abund_weight))%>%
	left_join(div.all,by='PLT_CN')
save(plots.div,file="data analysis2/plots.summary.rda")

## diveristy metric changes through time within invaded/uninvade plots
load("data analysis2/plots.summary.rda")
library(dplyr)
library(ggplot2);require(ggpmisc)
vars=c("PD.native","FAD.native","M_Ap.native","M_At.native")
index.labs <- c("PD", "FAD","MBL","MTD")
names(index.labs) <- vars
	
dat.plt=plots.div%>%filter(year>=2000&year<=2021)%>%select(pltID:disturb,all_of(vars))%>%
	tidyr::pivot_longer(cols=PD.native:M_At.native,names_to = "index", values_to = "diversity")%>%group_by(year,invaded,index)%>%
	summarize(mean=mean(diversity,na.rm=T),nplots=length(na.omit(diversity)),se=ifelse(nplots==0,NA,sd(diversity,na.rm=T)/sqrt(nplots)))
dat.plt$invaded=factor(dat.plt$invaded,levels=c(1,0))
dat.plt$index=factor(dat.plt$index,levels=vars)

pvals = do.call(rbind,lapply(levels(dat.plt$invaded), function(i) {
  do.call(rbind,lapply(levels(dat.plt$index), function(j) {   
      if(nrow(dat.plt[dat.plt$invaded==i & dat.plt$index==j, ]) > 1) {
        data.frame(invaded=i, index=j, 
                   p.val=coef(summary(lm(mean ~ year, data = dat.plt[dat.plt$invaded==i &dat.plt$index==j, ])))[2,4])
      }    
  }))
}))
# Keep only rows with p.val < 0.05
pvals = pvals[pvals$p.val < 0.05, ]
d.subset = dat.plt%>%filter(index%in%pvals$index & invaded %in% pvals$invaded) 

theme=theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),		
		strip.text = element_text(colour = 'black', size = rel(1.5)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2),
		legend.position= "right",
		legend.background = element_rect(fill = NA),
		legend.title=element_blank(),
		legend.text=element_text(face="bold",size=12))
ggplot(dat.plt, aes(x=year, y=mean,fill=invaded)) + 		
	geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se,color=invaded),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+		
	geom_point(shape=21,size=2,color="black",alpha=0.8) + 		
	labs(x="Year",y="Native tree diversity",fill="Plot type") +
	geom_smooth(aes(x = year, y = mean,color=invaded),data=d.subset,method="lm",linewidth=1,show.legend=FALSE,se =T,linetype=1,alpha=0.5)+
	#stat_poly_eq(aes(label =  paste(stat(p.value.label), "*\"\"",sep = ""),color=invaded),rr.digits = 2,p.digits = 2,label.x=1,label.y=c(1,0.92),size=5)+			
	scale_fill_manual(values=c("#F8766D","#619CFF"),labels=c("Invaded plots","Uninvaded plots"))+ 	
	scale_color_manual(values=c("#F8766D","#619CFF"),labels=c("Invaded plots","Uninvaded plots"))+ 		
	facet_wrap(~index, scales="free_y",ncol=2,labeller = labeller(index = index.labs))+			
	theme+ guides(fill = guide_legend(override.aes = list(size=5)))

mf=function(div,dat.plt,inv){
	dat=dat.plt%>%filter(index==div&invaded==inv)
	m=lm(mean~year,dat)
	re=c(div=div,inv=inv,F = round(summary(m)$f[1],2),p=round(summary(m)$coef[2,4],2))
	return(re)
}
rbind(do.call(rbind,lapply(c("PD.native","FAD.native","M_Ap.native","M_At.native"),mf,dat.plt,inv=1)), do.call(rbind,lapply(c("PD.native","FAD.native","M_Ap.native","M_At.native"),mf,dat.plt,inv=0)))


##################################################
## Regression on evaluating impact of non natives
##################################################

## Identify all FIA plots where new species (native or non-native) appeared between t1 and t2. 
## Also, the new species must still be present at t3 (otherwise, it is not a successful invasion, and we cannot quantify its effect).

# Response variables:
# Rate of change in diversity ("D", which could be richness, MAD, or MAP) from t2 to t3, quantified as (D3 - D2)/(t3 - t2), with positive values indicating an increase in diversity.
# Rate of change in plot biomass (Mg/ha/yr) from t2 to t3, quantified as (B3 - B2)/(t3 - t2). 
# I would exclude the invader (whether native or non-native) from the response variables (change in diversity and change in biomass), because I think this will lead to more easily interpretable results (especially for biomass). If desired, you could also consider an alternative diversity version, where you include the invader (this seems especially relevant for MAP, since non-native invaders might add substantial phylogenetic diversity).

# Explanatory variables (some of these are different from what I suggested before... it is not 100% clear what the best choices would be):
# EITHER a dummy variable indicating that the invader is native or non-native, OR a species random effect. I would not include both, because they will be confounded with each other if the within-group variance (within natives and within non-natives) is small compared to the between-group variance.
# diversity metric at t2 (I would only include the one that corresponds to the response variable)
# stand age at t2
# plot biomass at t2 (B2, Mg/ha)
# Percent of B2 that belongs to the invader, OR you could try the t2 biomass of the invader (Mg/ha). Whether you use percent or Mg/ha, the value will depend on the number of tallied individuals of the invading species, their TPHA values, and their individual biomasses. 
# ecoregion random effect
library(dplyr)
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/fia.abundance.rda")
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/traits.impute.rda")
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/phylo.imputate.Rdata")

# load("data analysis2/fia.abundance.rda")
# load("data analysis2/traits.impute.rda")
# load("data analysis2/phylo.imputate.Rdata")

status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=fia%>%filter(STATUSCD==1&LON>=-100&CONDPROP_UNADJ==1&(!degreeOfEstablishment%in%status[3]))
allplots$measure.time=allplots$MEASYEAR+allplots$MEASMON/12
sp.status=allplots%>%filter(!is.na(PREV_PLT_CN))%>%select(Accepted_name,degreeOfEstablishment)%>%distinct
sp.status$inv=ifelse(is.na(sp.status$degreeOfEstablishment),0,1)
#plots.cn=allplots%>%select(PLT_CN)

get.plt.reg=function(plt,allplots,sp.status,tree,trait_mat){
#re.all=c()
#for(plt in plots.cn$PLT_CN[which(plots.cn$PLT_CN==plt):nrow(plots.cn)]){
#dat=re.all%>%filter(Bcol.native>0&SR.t23>0)
# plt=dat$PLT_CN[1]
	require(dplyr);require(hillR);require(ape);require(FD)
	div.cal=function(abund_mat,tree,trait_mat){
	abund_mat[is.na(abund_mat)]=0
	comm =as.matrix(abund_mat)
	comm <- sweep(comm , 1, rowSums(abund_mat, na.rm = TRUE), "/") #abd(i)/sum(abd)
	comm <- comm[,comm > 0];
	df2.phy=abund_mat[,colnames(abund_mat)%in%tree$tip.label]
	if (length(df2.phy)<=1){
		M_p = M_Ap = PD = NA
	}else{
		tre.plt=keep.tip(tree, tip=names(df2.phy))
		PD=hill_phylo(df2.phy,tre.plt, q = 0)#equal to picante::pd(), do not account for abundance	
		dis2=cophenetic(tre.plt)#pairwise branch length
		df2.phy=as.matrix(df2.phy)
		M_Ap = sum(df2.phy %*% dis2 %*% matrix(df2.phy, ncol = 1))/(nrow(dis2)*(nrow(dis2)-1))	
	}
	if (length(comm)<=1){
		M_t = M_At = FAD =NA
	}else{	
		trait_mat.t=trait_mat[names(comm),]
		gowdis_corect=function(n) {
			n[which.max(n)]=ifelse(max(n,na.rm = T)==min(n,na.rm = T)&length(na.omit(n))>1,
			max(n,na.rm = T)+min(n,na.rm = T)/10000,max(n,na.rm = T));
			return(n)
		}
		trait_mat.t=apply(trait_mat.t,2,gowdis_corect)
		G.trait.mat=gowdis(trait_mat.t,ord = "podani")
		if(sum(G.trait.mat)>0){
			dis2 <- as.matrix(G.trait.mat)
			M_At = sum(comm %*% dis2 %*% matrix(comm, ncol = 1))/(nrow(dis2)*(nrow(dis2)-1))						
		}else{
			M_t = M_At = NA
		}
		FAD=sum(G.trait.mat)	
	}
	re=data.frame(M_Ap,M_At,PD,FAD)
	return(re)	
}

	t2=allplots%>%filter(PLT_CN%in%plt)
	t3=allplots%>%filter(PREV_PLT_CN%in%unique(t2$PLT_CN))
	if (nrow(t3)==0) return(NULL) # next; ## No plot mears or not used in t3
	new.inv.t2=t2%>%filter(is.na(PREV_TRE_CN))%>%distinct
	if (nrow(new.inv.t2)==0) return(NULL) #next; ## No new individuals in t2
	t1=allplots%>%filter(PLT_CN%in%unique(t2$PREV_PLT_CN))	
	if (nrow(t1)==0) return(NULL) #next; ## No plot mears or not used in t1
	
	new.sp.suvor=new.inv.t2%>%select(Accepted_name,CN,Biomass)%>%left_join(sp.status,by="Accepted_name")
	new.sp.suvor$suvor.t1=new.sp.suvor$Accepted_name%in%t1$Accepted_name #t2 new individuals of which species was present at t1
	new.sp.suvor$suvor.t3=new.sp.suvor$CN%in%t3$PREV_TRE_CN #t2 new individuals of which species survived to t3 
	
	suvor.t1=sum(new.sp.suvor$suvor.t1)/length(new.sp.suvor$suvor.t1)#proportion of t2 new individuals of which species was present at t1
	suvor.t1.biomass=sum(new.sp.suvor$Biomass[new.sp.suvor$suvor.t1])/sum(new.sp.suvor$Biomass)#proportion of t2 new individuals of which species was present at t1 (biomass)
	
	measure.time12=unique(t2$measure.time)-unique(t1$measure.time)
	measure.time23=unique(t3$measure.time)-unique(t2$measure.time)
	
	si.tot=unique(new.sp.suvor[new.sp.suvor$suvor.t1==F&new.sp.suvor$suvor.t3==T,"Accepted_name"])#successful invader (survival to t3) and the species was not present at t1.
	inv.tot=unique(new.sp.suvor[new.sp.suvor$suvor.t1==F,"Accepted_name"])#invaders whose the species was not present at t1, do not require survival to t3
	si.nn=unique(new.sp.suvor[new.sp.suvor$suvor.t1==F&new.sp.suvor$suvor.t3==T&new.sp.suvor$inv==1,"Accepted_name"])#successful non-native invader and the species was not present at t1.
	
	vars.t2=t2%>%group_by(PLT_CN,PREV_PLT_CN)%>%
		summarize(STDAGE=mean(STDAGE,na.rm=T),ECOSUBCD=unique(ECOSUBCD),
			Biomass.t2=sum(Biomass[(!CN%in%new.sp.suvor$CN)&is.na(degreeOfEstablishment)]),#biomass at t2, exclude t2 invaders and non-natives
			SR.t2=length(unique(Accepted_name[(!CN%in%new.sp.suvor$CN)&is.na(degreeOfEstablishment)])),#sp richness at t2,  exclude t2 invaders and non-natives
			BA.t2=sum(Abund_weight[(!CN%in%new.sp.suvor$CN)&is.na(degreeOfEstablishment)]),#basal area at t2,  exclude t2 invaders and non-natives
			nn.t2=sum(Biomass[!is.na(degreeOfEstablishment)]),#t2 biomass of non-natives, incl. both t2 invaders and t1 survivors
			si.tot.t12= sum(Biomass[Accepted_name%in%si.tot]), # t2 biomass of individuals that survived to t3, where the species is non-native or native and was not present at t1
			inv.tot.t12= sum(Biomass[Accepted_name%in%inv.tot]), # the same as si.tot.t12, except we do not require survival to t3			
			si.nn.t12.por=sum(Biomass[Accepted_name%in%si.nn])/si.tot.t12, # Proportion of successful invaders (si) biomass that belongs to non-native (nn) species
			
			SR.si.tot.t12= length(si.tot), # t2 SR of individuals that survived to t3, where the species is non-native or native and was not present at t1
			SR.inv.tot.t12= length(inv.tot), # the same as si.tot.t12, except we do not require survival to t3			
			SR.si.nn.t12.por=length(si.nn)/length(si.tot),
			SR.nn.t2=length(unique(Accepted_name[!is.na(degreeOfEstablishment)])),#t2 SR of non-natives, incl. both t2 invaders and t1 survivors			
			.groups="keep") # Proportion of successful invaders (si) richness that belongs to non-native (nn) species
	vars.t3=t3%>%group_by(PREV_PLT_CN)%>%
		summarize(Biomass.t3=sum(Biomass[(!PREV_TRE_CN%in%new.sp.suvor$CN)&is.na(degreeOfEstablishment)]),#biomass at t3, exclude t2 invaders and non-natives
			SR.t3=length(unique(Accepted_name[(!PREV_TRE_CN%in%new.sp.suvor$CN)&is.na(degreeOfEstablishment)])),#sp richness at t3,  exclude t2 invaders and non-natives
			BA.t3=sum(Abund_weight[(!PREV_TRE_CN%in%new.sp.suvor$CN)&is.na(degreeOfEstablishment)]),#basal area at t3,  exclude t2 invaders and non-natives			
			nn.t3=sum(Biomass[!is.na(degreeOfEstablishment)]),# t3 biomass of non-natives
			SR.nn.t3=length(unique(Accepted_name[!is.na(degreeOfEstablishment)])))# t3 SR of non-natives
	delta.t23=vars.t2%>%left_join(vars.t3,by=c("PLT_CN"="PREV_PLT_CN"))	%>%
		mutate(SR.t23=(SR.t3-SR.t2)/measure.time23, BA.t23=(BA.t3-BA.t2)/measure.time23, Biomass.t23=(Biomass.t3-Biomass.t2)/measure.time23,
			nn.t23=(nn.t3-nn.t2)/measure.time23, SR.nn.t23=(SR.nn.t3-SR.nn.t2)/measure.time23)	
	
	div.t2 <- t2%>%filter((!CN%in%new.sp.suvor$CN)&is.na(degreeOfEstablishment)) %>% group_by(PLT_CN,Accepted_name)  %>% summarize(abund_w = sum(Abund_weight),.groups ="keep")%>% 
		tidyr::pivot_wider(names_from = Accepted_name, values_from = abund_w) %>% tibble::column_to_rownames("PLT_CN")	%>%div.cal(tree,trait_mat)
	div.t3 <- t3%>%filter((!CN%in%new.sp.suvor$CN)&is.na(degreeOfEstablishment)) %>% group_by(PREV_PLT_CN,Accepted_name)  %>% summarize(abund_w = sum(Abund_weight),.groups ="keep")%>% 
		tidyr::pivot_wider(names_from = Accepted_name, values_from = abund_w) %>% tibble::column_to_rownames("PREV_PLT_CN")	%>%div.cal(tree,trait_mat)
	div.t23=(div.t3-div.t2)/measure.time23
	colnames(div.t23)=paste0(colnames(div.t23),".t23")
	re=data.frame(delta.t23,div.t2,div.t23,suvor.t1,suvor.t1.biomass,measure.time12,measure.time23)
	return(re)
	#re.all=rbind(re.all,re)
	#print(which(plots.cn$PLT_CN==plt))
#}
}

library(parallel)
no_cores <- 27#detectCores() - 1
mycl <- makePSOCKcluster(no_cores); 
plots.cn=allplots%>%filter(pltID.no.year%in%plots$pltID.no.year)%>%select(PLT_CN,PREV_PLT_CN)%>%distinct%>%na.omit
re.all=do.call(rbind,parLapply(cl=mycl,plots.cn$PLT_CN,get.plt.reg,allplots,sp.status,tree,trait_mat))
stopCluster(mycl)
save(re.all,file="/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.regression.rda")

library(dplyr)
# hist(re.all$measure.time12,breaks=20);abline(v=4,col="red");abline(v=6,col="red");abline(v=3,col="blue");abline(v=8,col="blue");
# windows()
# hist(re.all$measure.time23,breaks=20);abline(v=4,col="red");abline(v=6,col="red");abline(v=3,col="blue");abline(v=8,col="blue");
# plot(dat$si.tot.t12,dat$inv.tot.t12);abline(a=0,b=1,col="red")
# dat%>%filter(si.tot.t12==0&inv.tot.t12>10)%>%nrow
load("data analysis2/plots.regression2.rda")
library(lme4);
library(afex)#This will automatically add a p-value column to the output of the lmer for the fixed effects.
library(ggplot2)
re.all$Bcol.non_native=re.all$si.nn.t12.por*re.all$si.tot.t12
re.all$Bcol.native=re.all$si.tot.t12-re.all$Bcol.non_native
re.all$Bsur.non_native=ifelse(is.na(re.all$si.nn.t12.por),re.all$nn.t2,re.all$nn.t2-re.all$Bcol.non_native)

yvars=c("SR.t23","PD.t23","FAD.t23","M_Ap.t23","M_At.t23","Biomass.t23")
yvar.labs <- c("ΔSR", "ΔPD","ΔFAD","ΔMBL","ΔMTD","ΔBiomass")
names(yvar.labs) <- yvars

xvars=c("SR.t2","PD","FAD","M_Ap","M_Ap",NA)
stat.all=c()
for(i in 1:length(yvars)){
	yvar=yvars[i];xvar=xvars[i]
	dat=re.all%>%#filter(measure.time23>3&measure.time23<10&measure.time12>3&measure.time12<10)%>%
		filter(!is.na(Bcol.non_native))%>%filter(!(Bcol.non_native==0&Bcol.native==0&Bsur.non_native==0))
	#stat.n=dat%>%group_by(inv)%>%tally
	colnames(dat)[colnames(dat)==yvar]="yvar"
	if (!is.na(xvar)) {
		colnames(dat)[colnames(dat)==xvar]="Diversity.t2"
		# xvars.re=c("si.nn.t12.por","nn.t2","suvor.t1.biomass","si.tot.t12","Diversity.t2","Biomass.t2","STDAGE","measure.time12","measure.time23")
		# xvar.labs <- c("Colonizer.Non-native", "Non-native","Survivor","Colonizer","Diveristy","Biomass","StandAge","t12","t23")
		
		xvars.re=c("Bcol.non_native","Bsur.non_native","Bcol.native","Biomass.t2","Diversity.t2","STDAGE")
		xvar.labs <- c("Bcol.non_native","Bsur.non_native","Bcol.native","B2","Y2","StandAge")
		
		formu=as.formula(paste("yvar ~",paste(xvars.re,collapse="+"),"+(1|ECOSUBCD)",sep=" "))				
	}else{
		# xvars.re=c("si.nn.t12.por","nn.t2","suvor.t1.biomass","si.tot.t12","Biomass.t2","STDAGE","measure.time12","measure.time23")
		# xvar.labs <- c("Colonizer.Non-native", "Non-native","Survivor","Colonizer","Biomass","StandAge","t12","t23")
		xvars.re=c("Bcol.non_native","Bsur.non_native","Bcol.native","Biomass.t2","STDAGE")
		xvar.labs <- c("Bcol.non_native","Bsur.non_native","Bcol.native","B2","StandAge")
		
		formu=as.formula(paste("yvar ~",paste(xvars.re,collapse="+"),"+(1|ECOSUBCD)",sep=" "))		
	}
	dat[,xvars.re]=scale(dat[,xvars.re])
	M <- lmer(formu,data=dat)
	stat=as.data.frame(summary(M)$coef[-1,c(1,2,5)])
	colnames(stat)[2:3]=c("se","p")
	#rownames(stat)[rownames(stat)=="Diversity.t2"]=xvar
	stat$var=xvar.labs
	stat$var=factor(stat$var,levels=rev(xvar.labs))
	stat$type=ifelse((stat$Estimate+1.96*stat$se)<0,"negative",ifelse((stat$Estimate-1.96*stat$se)>0,"positive","insig"))
	stat$yvar=yvar
	stat.all=rbind(stat.all,stat)
}
stat.all$se=1.96*stat.all$se
stat.all2=stat.all%>%select(-type)%>%tidyr::pivot_wider(names_from = var, values_from = c(Estimate,se,p))%>%
	select(yvar,contains("Bcol.non_native"),contains("Bsur.non_native"),contains("Bcol.native"),contains("B2"),contains("Y2"),contains("StandAge"))	
# write.csv(stat.all2,"data analysis2/Fig.3.data.csv")
stat.all$yvar=factor(stat.all$yvar,levels=yvars)
ggplot() +geom_errorbar(data=stat.all,aes(x=Estimate, y=var,xmin=Estimate-se, xmax=Estimate+se,color=type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
	geom_point(data=stat.all,aes(x=Estimate, y=var,fill=type),size=2,shape=21,color="black",show.legend=F) + 		
	scale_fill_manual(values=c("insig"="darkgray", "positive"="#6D9EC1","negative"="#E46726"))+	
	scale_color_manual(values=c("insig"="darkgray", "positive"="#6D9EC1","negative"="#E46726"))+	
	labs(x="Estimates") +
	# geom_text(data = stat.n%>%filter(inv=="native"), mapping = aes(x=0 , y=-0.5, label = paste0("Plots (Native) = ",n)), hjust   = "center", vjust   = "inward",size=4)+
	# geom_text(data = stat.n%>%filter(inv=="alien"), mapping = aes(x=0 , y=-1.2, label = paste0("Plots (NonNative) = ",n)), hjust   = "center", vjust   = "inward",size=4)+	
	facet_wrap(~yvar, scales="free_x",labeller = labeller(yvar = yvar.labs))+		
	theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),
		axis.title.y=element_blank(),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))+
	geom_vline(xintercept=0,col="lightgray",size=0.8,linetype="longdash")

## new tree in t3
require(dplyr)
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/fia.abundance.rda")
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.regression.rda")
load("data analysis2/plots.regression.rda")
load("data analysis2/fia.abundance.rda")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=fia%>%filter(STATUSCD==1&LON>=-100&CONDPROP_UNADJ==1&(!degreeOfEstablishment%in%status[3]))
allplots$measure.time=allplots$MEASYEAR+allplots$MEASMON/12
sp.status=allplots%>%filter(!is.na(PREV_PLT_CN))%>%select(Accepted_name,degreeOfEstablishment)%>%distinct
sp.status$inv=ifelse(is.na(sp.status$degreeOfEstablishment),0,1)
re.all$Bcol.non_native=re.all$si.nn.t12.por*re.all$si.tot.t12
re.all$Bcol.native=re.all$si.tot.t12-re.all$Bcol.non_native
dat=re.all%>%filter(Bcol.native>0)%>%left_join(fia%>%select(PLT_CN,pltID.no.year,pltID)%>%distinct,by="PLT_CN")

get.nwtree=function(plt,allplots,sp.status){
	require(dplyr)
	t2=allplots%>%filter(PLT_CN%in%plt)
	t3=allplots%>%filter(PREV_PLT_CN%in%unique(t2$PLT_CN))
	if (nrow(t3)==0) return(NULL) # next; ## No plot mears or not used in t3
	new.inv.t2=t2%>%filter(is.na(PREV_TRE_CN))%>%distinct
	if (nrow(new.inv.t2)==0) return(NULL) #next; ## No new individuals in t2
	t1=allplots%>%filter(PLT_CN%in%unique(t2$PREV_PLT_CN))	
	if (nrow(t1)==0) return(NULL) #next; ## No plot mears or not used in t1
	new.tree.t3=t3%>%filter(is.na(PREV_TRE_CN)&(!Accepted_name%in%t2$Accepted_name))%>%distinct%>%left_join(sp.status,by="Accepted_name")%>%filter(inv==0)
	if (nrow(new.tree.t3)==0) return(NULL)
	sp.t3=new.tree.t3%>%filter(!Accepted_name%in%t1$Accepted_name)%>%select(Accepted_name)%>%nrow
	sp.t13=new.tree.t3%>%filter(Accepted_name%in%t1$Accepted_name)%>%select(Accepted_name)%>%nrow
	re=data.frame(PLT_CN=plt,sp.t3,sp.t13)
}
	
library(parallel)
no_cores <- detectCores() - 1
mycl <- makePSOCKcluster(no_cores); 
plots.cn=dat$PLT_CN
nw=do.call(rbind,parLapply(cl=mycl,plots.cn,get.nwtree,allplots,sp.status))
stopCluster(mycl)
dat=dat%>%left_join(nw,by="PLT_CN")

dat%>%filter(sp.t3>0)%>%summarize(n=length(sp.t3),nplots=length(unique(pltID)),mean=mean(sp.t3),se=sd(sp.t3)/sqrt(length(sp.t3)))
dat%>%filter(sp.t13>0)%>%summarize(n=length(sp.t13),nplots=length(unique(pltID)),mmean=mean(sp.t13),se=sd(sp.t13)/sqrt(length(sp.t13)))
dat2=dat%>%filter(sp.t13>0&sp.t3>0)
dat2$diff=dat2$sp.t3-dat2$sp.t13
sd(dat2$diff)/sqrt(length(dat2$diff))
######################################################
#are individuals that disappeared in t3 more likely to be small and rare?
######################################################
# Compare the dead trees in t3 among each case in their num. of individuals of the species, DIA, HT, BA, biomass and their differences with the biggese in the PLOT
library(dplyr)
# load("/blue/matthewthomas1/yunpeng.liu/data analysis2/fia.abundance.rda")
# load("/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.regression.rda")
# load("/blue/matthewthomas1/yunpeng.liu/data analysis2/traits.impute.rda")

load("data analysis2/plots.regression.rda")
load("data analysis2/fia.abundance.rda")
load("data analysis2/traits.impute.rda")
trait_mat[is.na(trait_mat)]=0
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=fia%>%filter(STATUSCD%in%c(1,2)&LON>=-100&CONDPROP_UNADJ==1&(!degreeOfEstablishment%in%status[3]))#STATUSCD==1

#allplots%>%group_by(STATUSCD)%>%tally
#plt=re.all$PLT_CN[12]
get.dead.rare=function(plt,allplots,trait_mat){
	require(dplyr);require(FD)
	t2=allplots%>%filter(PLT_CN%in%plt)
	t3=allplots%>%filter(PREV_PLT_CN%in%unique(t2$PLT_CN))
	dead.t3=t3%>%filter(STATUSCD==2&(!is.na(PREV_TRE_CN))&is.na(degreeOfEstablishment)) #native tree dead in t3, including 1) tree live in t2 while dead in t3; 2)tree dead in t2 and still there in t3
	live.t3=t3%>%filter(STATUSCD==1&(!is.na(PREV_TRE_CN))&is.na(degreeOfEstablishment)) # t2 native trees survived in t3
	if (nrow(dead.t3)==0|nrow(live.t3)==0) return(NULL)
	dead.t3=dead.t3%>%filter(!Accepted_name%in%live.t3$Accepted_name)#dead tree that belong to species dead in t3
	if (nrow(dead.t3)==0) return(NULL)
	dead.t2=t2%>%filter(CN%in%dead.t3$PREV_TRE_CN&STATUSCD==1)#native tree live in t2 while dead in t3
	if (nrow(dead.t2)==0) return(NULL)
		
	# ## caculate Abundance weighted functional distinctiveness
	# dead.sp=unique(dead.t2$Accepted_name)
	# native.sp=t2%>%filter(is.na(degreeOfEstablishment))%>%select(Accepted_name)%>%distinct
	# dead.sp=dead.sp[dead.sp%in%native.sp$Accepted_name]
	# SR.loss=length(dead.sp)	
	# if(SR.loss==0) return(NULL) #in plt=="102139660010661", dead.sp="Morus alba" in t2 while it is "Morus rubra" in t3
	# abd=t2%>%group_by(Accepted_name)%>%tally
	# trait_mat.t=trait_mat[native.sp$Accepted_name,]
	# gowdis_corect=function(n) {
		# n[which.max(n)]=ifelse(max(n,na.rm = T)==min(n,na.rm = T)&length(na.omit(n))>1,
		# max(n,na.rm = T)+min(n,na.rm = T)/10000,max(n,na.rm = T));
		# return(n)
	# }
	# trait_mat.t=apply(trait_mat.t,2,gowdis_corect)
	# G.trait.mat=gowdis(trait_mat.t,ord = "podani")%>%as.matrix
	# mta.cal=function(i.sp,G.trait.mat,abd,sp.all){
		# j.sp=sp.all[!sp.all%in%i.sp]
		# if (length(j.sp)==0) {
			# return(NA)
		# }else{
			# dij=as.matrix(G.trait.mat[j.sp,i.sp]);colnames(dij)=i.sp
			# nj=abd$n[abd$Accepted_name%in%j.sp];names(nj)=j.sp
			# MTA=sum(dij*nj)/sum(nj)
			# return(MTA)
		# }		
	# }
	# mta.dead=do.call(mean,lapply(dead.sp,mta.cal,G.trait.mat,abd,sp.all=native.sp$Accepted_name))
		
	dorminant=t2%>%filter(STATUSCD==1)%>%group_by(Accepted_name)%>%summarise(Abd=sum(TPA_UNADJ),Biomass=sum(Biomass,na.rm = T))
	dorminant$Abd=dorminant$Abd/sum(dorminant$Abd);dorminant$Biomass=dorminant$Biomass/sum(dorminant$Biomass,na.rm = T)
	rare.Abd=dorminant[dorminant$Abd<0.2,"Accepted_name"] #rare species in a plot is defined as species whose total BA < 20% of overall BA in the plot
	rare.biomass=dorminant[dorminant$Biomass<0.2,"Accepted_name"] #rare species in a plot is defined as species whose total Biomass < 20% of overall Biomass in the plot
	#max.cond=t2%>%filter(STATUSCD==1)%>%group_by(PLT_CN)%>%summarise(BA.max=max(Abund_weight),Biomass.max=max(Biomass,na.rm=T),DIA.max=max(DIA,na.rm=T),HT.max=max(HT,na.rm=T))
	#dead.t2.stat=dead.t2%>%group_by(PLT_CN,Accepted_name,CN)%>%summarise(BA=Abund_weight/max.cond$BA.max,Biomass=Biomass/max.cond$Biomass.max,DIA=DIA/max.cond$DIA.max,HT=HT/max.cond$HT.max,.groups="keep")
	
	rare.individual.por.ba=length(dead.t2$Accepted_name[dead.t2$Accepted_name%in%rare.ba$Accepted_name])/length(dead.t2$Accepted_name)# proportion of dead individuals in t3 that are rare species (BA)
	rare.sp.por.ba=length(unique(dead.t2$Accepted_name[dead.t2$Accepted_name%in%rare.ba$Accepted_name]))/length(unique(dead.t2$Accepted_name))# proportion of rare species to all the species dead in t3 (BA)

	# rare.individual.por=length(dead.t2$Accepted_name[dead.t2$Accepted_name%in%rare.biomass$Accepted_name])/length(dead.t2$Accepted_name)# proportion of dead individuals in t3 that are rare species (Biomass)
	# rare.sp.por=length(unique(dead.t2$Accepted_name[dead.t2$Accepted_name%in%rare.biomass$Accepted_name]))/length(unique(dead.t2$Accepted_name))# proportion of rare species to all the species dead in t3 (Biomass)

	# rare.Biomass.por.ba=sum(dead.t2$Biomass[dead.t2$Accepted_name%in%rare.ba$Accepted_name],na.rm=T)/sum(dead.t2$Biomass,na.rm=T)# proportion of the Biomass of dead individuals in t3 that are rare species (BA)
	# rare.Biomass.por=sum(dead.t2$Biomass[dead.t2$Accepted_name%in%rare.biomass$Accepted_name],na.rm=T)/sum(dead.t2$Biomass,na.rm=T)# proportion of the Biomass of rare species to all the species dead in t3 (Biomass)

	# rare.BA.por.ba=sum(dead.t2$Abund_weight[dead.t2$Accepted_name%in%rare.ba$Accepted_name])/sum(dead.t2$Abund_weight)# proportion of the BA of dead individuals in t3 that are rare species (BA)
	# rare.BA.por=sum(dead.t2$Abund_weight[dead.t2$Accepted_name%in%rare.biomass$Accepted_name])/sum(dead.t2$Abund_weight)# proportion of the BA of rare species to all the species dead in t3 (Biomass)
	re=data.frame(PLT_CN=plt,rare.individual.por.ba,rare.sp.por.ba)
	# re=data.frame(PLT_CN=plt,SR.loss,mta.dead,
		# rare.individual.por.ba,rare.sp.por.ba,rare.Biomass.por.ba,rare.BA.por.ba,rare.individual.por,rare.sp.por,rare.Biomass.por,rare.BA.por)
	return(re)
}
library(parallel)
no_cores <- 37#detectCores() - 1
mycl <- makePSOCKcluster(no_cores); 
re.all.dead=do.call(rbind,parLapply(cl=mycl,re.all$PLT_CN,get.dead.rare,allplots,trait_mat))
stopCluster(mycl)
save(re.all.dead,file="/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.dead.trees.V3.rda")

# re.all.dead=c()	
# for (i in re.all$PLT_CN){
	# tmp=get.dead.rare(i,allplots,trait_mat)
	# if(!is.null(tmp)) re.all.dead=rbind(re.all.dead,tmp)
	# print(which(re.all$PLT_CN==i))
# }

######################################################
## Is the colonizer displacing an ecologically similar species?
######################################################
# compute distinctiveness for all species in each community (eq. 4 in REF: Conservation prioritization based on trait based metrics illustrated with global parrot distributions)
require(dplyr)
# load("/blue/matthewthomas1/yunpeng.liu/data analysis2/fia.abundance.rda")
# load("/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.regression.rda")
# load("/blue/matthewthomas1/yunpeng.liu/data analysis2/traits.impute.rda")
# load("/blue/matthewthomas1/yunpeng.liu/data analysis2/phylo.imputate.Rdata")

load("data analysis2/traits.impute.rda")
load("data analysis2/fia.abundance.rda")
load("data analysis2/plots.regression.rda")
load("data analysis2/phylo.imputate.Rdata")

plt.all=re.all%>%filter(nn.t2!=0)
trait_mat[is.na(trait_mat)]=0
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=fia%>%filter(STATUSCD%in%c(1,2)&LON>=-100&CONDPROP_UNADJ==1&(!degreeOfEstablishment%in%status[3]))#STATUSCD==1
sp.status=fia%>%select(Accepted_name,degreeOfEstablishment)%>%distinct
sp.status$inv=ifelse(is.na(sp.status$degreeOfEstablishment),0,1)
#continetal scale rareness
native=allplots%>%filter(STATUSCD==1)%>%group_by(Accepted_name)%>%summarise(Abd=sum(TPA_UNADJ),Biomass=sum(Biomass,na.rm = T))%>%as.data.frame
native$Abd=native$Abd/sum(native$Abd);native$Biomass=native$Biomass/sum(native$Biomass,na.rm = T)
rare.all=list(rare1.1.all=native[native$Abd<0.01,"Accepted_name"],
	rare1.2.all=native[native$Abd<0.05&native$Abd>=0.01,"Accepted_name"],
	rare1.3.all=native[native$Abd<=0.1&native$Abd>=0.05,"Accepted_name"],
	rare1.4.all=native[native$Abd>0.1,"Accepted_name"],
	rare2.1.all=native[native$Biomass<0.01,"Accepted_name"] ,
	rare2.2.all=native[native$Biomass<0.05&native$Biomass>=0.01,"Accepted_name"],
	rare2.3.all=native[native$Biomass<=0.1&native$Biomass>=0.05,"Accepted_name"],
	rare2.4.all=native[native$Biomass>0.1,"Accepted_name"])
																						
get.mta2=function(plt,trait_mat=NULL,tree=NULL,allplots,sp.status,analysis="FDis",rare.all=NULL){
	#plt=plt.all$PLT_CN[1]
	require(dplyr)
	t2=allplots%>%filter(PLT_CN%in%plt)
	t1=allplots%>%filter(PLT_CN%in%unique(t2$PREV_PLT_CN))	
	if (nrow(t1)==0) return(NULL) #next; ## No plot mears or not used in t1	
	t3=allplots%>%filter(PREV_PLT_CN%in%unique(t2$PLT_CN))
	dead.t3=t3%>%filter(STATUSCD==2&(!is.na(PREV_TRE_CN))) #tree dead in t3, including 1) tree live in t2 while dead in t3; 2)tree dead in t2 and still there in t3
	if (nrow(dead.t3)==0) return(NULL)
	dead.t2=t2%>%filter(CN%in%dead.t3$PREV_TRE_CN&STATUSCD==1&is.na(degreeOfEstablishment))#&(!is.na(PREV_TRE_CN))#native tree live in t1,t2 while dead in t3
	if (nrow(dead.t2)==0) return(NULL)	
	new.inv.t2=t2%>%filter(is.na(PREV_TRE_CN))%>%distinct
	if (nrow(new.inv.t2)==0) return(NULL) #next; ## No new individuals in t2
	sp.all=unique(t2$Accepted_name)
	if (length(sp.all)<2)	return(NULL)
	
	new.sp.suvor=new.inv.t2%>%select(Accepted_name,CN,Biomass)%>%left_join(sp.status,by="Accepted_name")
	new.sp.suvor$suvor.t1=new.sp.suvor$Accepted_name%in%t1$Accepted_name #t2 new individuals of which species was present at t1
	#new.sp.suvor$suvor.t3=new.sp.suvor$CN%in%t3[t3$STATUSCD==1,"PREV_TRE_CN"] #t2 new individuals of which species survived to t3 
		
	#si.native=unique(new.sp.suvor[new.sp.suvor$suvor.t1==F&new.sp.suvor$suvor.t3==T&new.sp.suvor$inv==0,"Accepted_name"])#successful non-native invader and the species was not present at t1.
	nn=sp.all[sp.all%in%sp.status[sp.status$inv==1,"Accepted_name"]]
	si.nn=unique(new.sp.suvor[new.sp.suvor$suvor.t1==F&new.sp.suvor$inv==1,"Accepted_name"])#successful non-native invader and the species was not present at t1.
	dead.sp=unique(dead.t2$Accepted_name)
	native.survior=t2%>%filter(STATUSCD==1&(!CN%in%dead.t2$CN)&(!Accepted_name%in%sp.status[sp.status$inv==1,"Accepted_name"])&(!is.na(PREV_TRE_CN)))%>%
		as.data.frame%>%.[,"Accepted_name"]%>%unique#native species survived to t3
	dead.sp=dead.sp[!dead.sp%in%native.survior]
	if (analysis=="PhyloDis"){
		require(ape)
		sp.all=sp.all[sp.all%in%tree$tip]
		dead.sp=dead.sp[dead.sp%in%tree$tip]
		native.survior=native.survior[native.survior%in%tree$tip]
		if (length(sp.all)<2)	return(NULL)
		tre.plt=keep.tip(tree, tip=sp.all)
		phydis=stats::cophenetic(tre.plt)#pairwise branch length
		abd0=t2%>%group_by(Accepted_name)%>%summarise(n=sum(Abund_weight))%>%as.data.frame
		abd=abd0$n;names(abd)=abd0$Accepted_name
		mpa.cal=function(i.sp,phydis,abd,j.sp=NULL,sp.all=NULL){
			if (is.null(j.sp)) j.sp=sp.all[!sp.all%in%i.sp] 
			#if (is.null(j.sp)) j.sp=sp.all else j.sp=c(j.sp,i.sp)
			if (length(j.sp)==0) {
				return(NA)
			}else{
				dij=as.matrix(phydis[j.sp,i.sp]);colnames(dij)=i.sp
				nj=abd[j.sp];names(nj)=j.sp
				MPA=sum(dij*nj)/sum(nj)
				return(MPA)
			}
		}
		re=data.frame(PLT_CN=plt,
			mpa.dead=ifelse(length(dead.sp)==0,NA,do.call(c,lapply(dead.sp,mpa.cal,phydis,abd,j.sp=NULL,sp.all=sp.all))%>%
				stats::weighted.mean(.,abd[dead.sp],na.rm=T)),
			mpa.dead2=ifelse(length(dead.sp)==0,NA,do.call(c,lapply(dead.sp,mpa.cal,phydis,abd,j.sp=native.survior,sp.all=sp.all))%>%
				stats::weighted.mean(.,abd[dead.sp],na.rm=T)),
			mpa.native.survior=ifelse(length(native.survior)==0,NA,do.call(c,lapply(native.survior,mpa.cal,phydis,abd,j.sp=NULL,sp.all=sp.all))%>%
				stats::weighted.mean(.,abd[native.survior],na.rm=T)),
			mpa.native.survior2=ifelse(length(native.survior)==0,NA,do.call(c,lapply(native.survior,mpa.cal,phydis,abd,j.sp=native.survior,sp.all=sp.all))%>%
				stats::weighted.mean(.,abd[native.survior],na.rm=T))		
		)
		return(re)	
		
	}
	if (analysis=="FDis"){
		require(FD)
		## caculate Abundance weighted functional distinctiveness
		abd0=t2%>%group_by(Accepted_name)%>%summarise(n=sum(Abund_weight))%>%as.data.frame
		abd=abd0$n;names(abd)=abd0$Accepted_name
		trait_mat.t=trait_mat[sp.all,]
		gowdis_corect=function(n) {
			n[which.max(n)]=ifelse(max(n,na.rm = T)==min(n,na.rm = T)&length(na.omit(n))>1,
			max(n,na.rm = T)+min(n,na.rm = T)/10000,max(n,na.rm = T));
			return(n)
		}
		trait_mat.t=apply(trait_mat.t,2,gowdis_corect)
		G.trait.mat=gowdis(trait_mat.t,ord = "podani")%>%as.matrix
		mta.cal=function(i.sp,G.trait.mat,abd,j.sp=NULL,sp.all=NULL){
			if (is.null(j.sp)) j.sp=sp.all[!sp.all%in%i.sp] 
			#if (is.null(j.sp)) j.sp=sp.all else j.sp=c(j.sp,i.sp)
			if (length(j.sp)==0) {
				return(NA)
			}else{
				dij=as.matrix(G.trait.mat[j.sp,i.sp]);colnames(dij)=i.sp
				nj=abd$n[abd$Accepted_name%in%j.sp];names(nj)=j.sp
				MTA=sum(dij*nj)/sum(nj)
				return(MTA)
			}
		}
		re=data.frame(PLT_CN=plt,
			mta.si.nn=ifelse(length(si.nn)==0,NA,do.call(c,lapply(si.nn,mta.cal,G.trait.mat,abd,j.sp=native.survior,sp.all=sp.all))%>%
				stats::weighted.mean(.,abd[si.nn]$n,na.rm=T)),
			mta.nn=ifelse(length(nn)==0,NA,do.call(c,lapply(nn,mta.cal,G.trait.mat,abd,j.sp=native.survior,sp.all=sp.all))%>%
				stats::weighted.mean(.,abd[nn]$n,na.rm=T)),
			mta.dead=ifelse(length(dead.sp)==0,NA,do.call(c,lapply(dead.sp,mta.cal,G.trait.mat,abd,j.sp=native.survior,sp.all=sp.all))%>%
				stats::weighted.mean(.,abd[dead.sp]$n,na.rm=T)),
			mta.native.survior=ifelse(length(native.survior)==0,NA,do.call(c,lapply(native.survior,mta.cal,G.trait.mat,abd,j.sp=native.survior,sp.all=sp.all))%>%
				stats::weighted.mean(.,abd[native.survior]$n,na.rm=T)),	
			mta.si.nn.dead=ifelse(length(si.nn)==0|length(dead.sp)==0,NA,do.call(c,lapply(si.nn,mta.cal,G.trait.mat,abd,j.sp=dead.sp,sp.all=sp.all))%>%
				stats::weighted.mean(.,abd[si.nn]$n,na.rm=T)),
			mta.nn.dead=ifelse(length(nn)==0|length(dead.sp)==0,NA,do.call(c,lapply(nn,mta.cal,G.trait.mat,abd,j.sp=dead.sp,sp.all=sp.all))%>%
				stats::weighted.mean(.,abd[nn]$n,na.rm=T))
		)
		return(re)		
	}
	if (analysis=="rare"){
		dead.sp=dead.sp[!dead.sp%in%native.survior]
		if (length(native.survior)>0|length(dead.sp)>0){
			dorminant=t2%>%filter(STATUSCD==1)%>%group_by(Accepted_name)%>%summarise(Abd=sum(TPA_UNADJ),Biomass=sum(Biomass,na.rm = T))%>%as.data.frame
			dorminant$Abd=dorminant$Abd/sum(dorminant$Abd);dorminant$Biomass=dorminant$Biomass/sum(dorminant$Biomass,na.rm = T)
			rare1.1=dorminant[dorminant$Abd<0.01,"Accepted_name"];
			cat1.1.survior=length(native.survior[native.survior%in%rare1.1])
			cat1.1.dead=length(dead.sp[dead.sp%in%rare1.1])
			all.cat1.1.survior=length(native.survior[native.survior%in%rare.all[[1]]])
			all.cat1.1.dead=length(dead.sp[dead.sp%in%rare.all[[1]]])
			
			rare1.2=dorminant[dorminant$Abd<0.05&dorminant$Abd>=0.01,"Accepted_name"]
			cat1.2.survior=length(native.survior[native.survior%in%rare1.2])
			cat1.2.dead=length(dead.sp[dead.sp%in%rare1.2])
			all.cat1.2.survior=length(native.survior[native.survior%in%rare.all[[2]]])
			all.cat1.2.dead=length(dead.sp[dead.sp%in%rare.all[[2]]])
			
			rare1.3=dorminant[dorminant$Abd<=0.1&dorminant$Abd>=0.05,"Accepted_name"]
			cat1.3.survior=length(native.survior[native.survior%in%rare1.3])
			cat1.3.dead=length(dead.sp[dead.sp%in%rare1.3])
			all.cat1.3.survior=length(native.survior[native.survior%in%rare.all[[3]]])
			all.cat1.3.dead=length(dead.sp[dead.sp%in%rare.all[[3]]])
			
			rare1.4=dorminant[dorminant$Abd>0.1,"Accepted_name"]
			cat1.4.survior=length(native.survior[native.survior%in%rare1.4])
			cat1.4.dead=length(dead.sp[dead.sp%in%rare1.4])
			all.cat1.4.survior=length(native.survior[native.survior%in%rare.all[[4]]])
			all.cat1.4.dead=length(dead.sp[dead.sp%in%rare.all[[4]]])
			
			rare2.1=dorminant[dorminant$Biomass<0.01,"Accepted_name"] 
			cat2.1.survior=length(native.survior[native.survior%in%rare2.1])
			cat2.1.dead=length(dead.sp[dead.sp%in%rare2.1])
			all.cat2.1.survior=length(native.survior[native.survior%in%rare.all[[5]]])
			all.cat2.1.dead=length(dead.sp[dead.sp%in%rare.all[[5]]])
			
			rare2.2=dorminant[dorminant$Biomass<0.05&dorminant$Biomass>=0.01,"Accepted_name"]
			cat2.2.survior=length(native.survior[native.survior%in%rare2.2])
			cat2.2.dead=length(dead.sp[dead.sp%in%rare2.2])
			all.cat2.2.survior=length(native.survior[native.survior%in%rare.all[[6]]])
			all.cat2.2.dead=length(dead.sp[dead.sp%in%rare.all[[6]]])
						
			rare2.3=dorminant[dorminant$Biomass<=0.1&dorminant$Biomass>=0.05,"Accepted_name"]
			cat2.3.survior=length(native.survior[native.survior%in%rare2.3])
			cat2.3.dead=length(dead.sp[dead.sp%in%rare2.3])
			all.cat2.3.survior=length(native.survior[native.survior%in%rare.all[[7]]])
			all.cat2.3.dead=length(dead.sp[dead.sp%in%rare.all[[7]]])
			
			rare2.4=dorminant[dorminant$Biomass>0.1,"Accepted_name"]
			cat2.4.survior=length(native.survior[native.survior%in%rare2.4])
			cat2.4.dead=length(dead.sp[dead.sp%in%rare2.4])
			all.cat2.4.survior=length(native.survior[native.survior%in%rare.all[[8]]])
			all.cat2.4.dead=length(dead.sp[dead.sp%in%rare.all[[8]]])
			
			re=rbind(data.frame(type1="Density <1%",type2="Native surviors",persent=cat1.1.survior,persent.all=all.cat1.1.survior),
				data.frame(type1="Density <1%",type2="Extinct natives",persent=cat1.1.dead,persent.all=all.cat1.1.dead),
				data.frame(type1="Density <5%",type2="Native surviors",persent=cat1.2.survior,persent.all=all.cat1.2.survior),
				data.frame(type1="Density <5%",type2="Extinct natives",persent=cat1.2.dead,persent.all=all.cat1.2.dead),
				data.frame(type1="Density <10%",type2="Native surviors",persent=cat1.3.survior,persent.all=all.cat1.3.survior),
				data.frame(type1="Density <10%",type2="Extinct natives",persent=cat1.3.dead,persent.all=all.cat1.3.dead),
				data.frame(type1="Density >=10%",type2="Native surviors",persent=cat1.4.survior,persent.all=all.cat1.4.survior),
				data.frame(type1="Density >=10%",type2="Extinct natives",persent=cat1.4.dead,persent.all=all.cat1.4.dead),
				
				data.frame(type1="Biomass <1%",type2="Native surviors",persent=cat2.1.survior,persent.all=all.cat2.1.survior),
				data.frame(type1="Biomass <1%",type2="Extinct natives",persent=cat2.1.dead,persent.all=all.cat2.1.dead),
				data.frame(type1="Biomass <5%",type2="Native surviors",persent=cat2.2.survior,persent.all=all.cat2.2.survior),
				data.frame(type1="Biomass <5%",type2="Extinct natives",persent=cat2.2.dead,persent.all=all.cat2.2.dead),
				data.frame(type1="Biomass <10%",type2="Native surviors",persent=cat2.3.survior,persent.all=all.cat2.3.survior),
				data.frame(type1="Biomass <10%",type2="Extinct natives",persent=cat2.3.dead,persent.all=all.cat2.3.dead),
				data.frame(type1="Biomass >=10%",type2="Native surviors",persent=cat2.4.survior,persent.all=all.cat2.4.survior),
				data.frame(type1="Biomass >=10%",type2="Extinct natives",persent=cat2.4.dead,persent.all=all.cat2.4.dead)			
			)
			re$PLT_CN=plt
			re$sur.sp=length(native.survior);re$dead.sp=length(dead.sp)
			return(re)
		}else{
			return(NULL)
		}	
	}
}
# MTA=c()
# for (i in plt.all$PLT_CN){
# tmp=get.mta2(i,trait_mat,allplots,sp.status)
# get.mta2(plt,trait_mat=NULL,tree=tree,allplots,sp.status,analysis="PhyloDis")
# MTA=rbind(MTA,tmp)
# }
library(parallel)
no_cores <- detectCores() - 1
mycl <- makePSOCKcluster(no_cores);
MTA=do.call(rbind,parLapply(cl=mycl,plt.all$PLT_CN,get.mta2,trait_mat=trait_mat,allplots,sp.status))
MPA=do.call(rbind,parLapply(cl=mycl,plt.all$PLT_CN,get.mta2,trait_mat=NULL,tree=tree,allplots,sp.status,analysis="PhyloDis"))
rare=do.call(rbind,parLapply(cl=mycl,plt.all$PLT_CN,get.mta2,trait_mat=NULL,tree=NULL,allplots,sp.status,"rare",rare.all=rare.all))
stopCluster(mycl)
#save(MTA,file="/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.MTA.rda")
#save(MTA,file="data analysis2/plots.MTA2.rda")
#save(rare,file="/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.rare2.rda")#plot in re.all, 34674 plot mears (22873 plots) 
#save(rare,file="/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.rare.rda")#plot in plt.all, 1024 plot mears
save(MPA,file="/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.MPA.rda")

#plot rareness group
load("data analysis2/plots.rare.rda")
#load("data analysis2/plots.rare2.rda")
library(ggpubr)
all.sp=rare%>%select(PLT_CN,dead.sp,sur.sp)%>%distinct%>%summarize(dead=sum(dead.sp),sur=sum(sur.sp))
rare.stat=rare%>%na.omit%>%group_by(type1,type2)%>%summarize(value=sum(persent),value.all=sum(persent.all))
rare.stat$Type=ifelse(grepl("Biomass", rare.stat$type1),"Biomass (Mg/ha)","Density (trees/ha)")
mf=function(st) stringr::str_split(st," ")[[1]][2]
rare.stat$type1=do.call(c,lapply(rare.stat$type1,mf))
rare.stat$type1=factor(rare.stat$type1,levels=c("<1%","<5%","<10%",">=10%"))
rare.stat$value=as.numeric(rare.stat$value)
rare.stat$value.all=as.numeric(rare.stat$value.all)
rare.stat[rare.stat$type2=="Extinct natives","value"]=rare.stat[rare.stat$type2=="Extinct natives","value"]/all.sp$dead*100
rare.stat[rare.stat$type2=="Native surviors","value"]=rare.stat[rare.stat$type2=="Native surviors","value"]/all.sp$sur*100
rare.stat[rare.stat$type2=="Extinct natives","value.all"]=rare.stat[rare.stat$type2=="Extinct natives","value.all"]/all.sp$dead*100
rare.stat[rare.stat$type2=="Native surviors","value.all"]=rare.stat[rare.stat$type2=="Native surviors","value.all"]/all.sp$sur*100
p1=ggplot(rare.stat,aes(x=type1, y=value,fill=type2)) +
    geom_bar(position="dodge", stat="identity",show.legend=T)+
	#coord_cartesian(ylim=c(0.386,0.404))+
   # geom_errorbar(aes(x=type1, ymin=value*100-se*196, ymax=value*100+se*196,colour=type2), width=0.4,  alpha=0.9, size=1.3,position=position_dodge(.9),show.legend=F)+
	labs(x="Percent of total biomass or density\n within each plot",y="Frequency (%)",fill="Type of Species")+
	scale_x_discrete(labels=c("<1%" = "<1%","<5%"="1-5%","<10%"="5-10%",">=10%"=">=10%"))+	
	facet_wrap(~Type, scales="free_x")+
	ylim(0,62)+
	theme_classic()+theme(axis.text = element_text(size=12,color='black'),
		axis.text.x =element_text(angle=90),
		axis.title = element_text(size=12,color='black'),
		#axis.title.x =element_blank(),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2),
		legend.background=element_rect(fill='transparent'),
		legend.text=element_text(size=12),
		legend.title=element_text(face="bold",size=12),
		legend.position=c(0.7,0.85))
p2=ggplot(rare.stat,aes(x=type1, y=value.all,fill=type2)) +
    geom_bar(position="dodge", stat="identity",show.legend=F)+
	#coord_cartesian(ylim=c(0.386,0.404))+
   # geom_errorbar(aes(x=type1, ymin=value*100-se*196, ymax=value*100+se*196,colour=type2), width=0.4,  alpha=0.9, size=1.3,position=position_dodge(.9),show.legend=F)+
	labs(x="Percent of total biomass or density\n across the eastern USA",y="Frequency (%)",fill="Type of Species")+
	scale_x_discrete(labels=c("<1%" = "<1%","<5%"="1-5%","<10%"="5-10%",">=10%"=">=10%"))+	
	facet_wrap(~Type, scales="free_x")+
	ylim(0,62)+
	theme_classic()+theme(axis.text = element_text(size=12,color='black'),
		axis.text.x =element_text(angle=90),
		axis.title = element_text(size=12,color='black'),
		axis.text.y=element_blank(),
		axis.title.y=element_blank(),
		#strip.text = element_blank(),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2)
		)
ggarrange(p1,p2,ncol=2,align="h",labels="auto",hjust=0,vjust=1.5)
	
# plot the Difference
# ref:http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
load("data analysis2/plots.MTA2.rda")
library(ggpubr);library(dplyr)
dat=MTA %>% select(PLT_CN,mta.dead,mta.native.survior) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")#%>%filter(PLT_CN%in%plt$PLT_CN)
dat1=MTA %>% select(PLT_CN,mta.native.survior,mta.si.nn) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")
dat2=MTA %>% select(PLT_CN,mta.si.nn.dead,mta.si.nn) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")
dat3=MTA %>% select(PLT_CN,mta.si.nn.dead,mta.native.survior) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")

#dat3=MTA %>% select(PLT_CN,mta.native.survior,mta.nn,mta.nn.dead) %>%na.omit%>%tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")
t.test(MTA ~ Type,data= dat,paired = T)
t.test(MTA ~ Type,data= dat1,paired = T)
t.test(MTA ~ Type,data= dat2,paired = T)
t.test(MTA ~ Type,data= dat3,paired = T)

dat1=MTA %>% select(PLT_CN,mta.native.survior,mta.nn) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")
dat2=MTA %>% select(PLT_CN,mta.nn.dead,mta.nn) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")
dat3=MTA %>% select(PLT_CN,mta.nn.dead,mta.native.survior) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")


#PLOT
library(ggplot2)
plt2=rbind(dat,dat1,dat2,dat3)%>%group_by(Type)%>%summarise(value=mean(MTA),se=sd(MTA)/sqrt(length(MTA)));
plt2$Type=factor(plt2$Type,levels=c("mta.native.survior","mta.dead","mta.si.nn","mta.si.nn.dead"))
plt2$Type=factor(plt2$Type,levels=c("mta.native.survior","mta.dead","mta.nn","mta.nn.dead"))
ggplot(plt2) +
    #geom_point(aes(x=Type, y=value),shape=21,fill="#F8766D",color="black",size=3) +
	geom_bar( aes(x=Type, y=value), stat="identity", fill="gray", alpha=0.7)+
	coord_cartesian(ylim=c(0.3,0.48))+
    geom_errorbar(aes(x=Type, ymin=value-1.96*se, ymax=value+1.96*se), width=0.4, colour="black", alpha=0.9, size=1.3)+
	labs(x="Type of species",y="Mean functional distance")+
	#scale_x_discrete(labels=c("mta.dead" = "Extinct natives vs.\nNative surviors","mta.native.survior" = "Native surviors vs.\nThemselves", "mta.si.nn" = "Non-native colonizers vs.\nNative surviors","mta.si.nn.dead"="Non-native colonizers vs.\nExtinct natives"))+
	scale_x_discrete(labels=c("mta.dead" = "Extinct natives vs.\nNative surviors","mta.native.survior" = "Native surviors vs.\nThemselves", "mta.nn" = "Non-natives vs.\nNative surviors","mta.nn.dead"="Non-natives vs.\nExtinct natives"))+	
		theme_classic()+theme(axis.title.x=element_blank(),
		axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=12,color='black'))+
	annotate("text", x = 1, y = 0.46, label = "a", size = 8, color = "black")+
	annotate("text", x = 2, y = 0.46, label = "a", size = 8, color = "black")+
	annotate("text", x = 3, y = 0.46, label = "b", size = 8, color = "black")+
	annotate("text", x = 4, y = 0.46, label = "a", size = 8, color = "black")

load("data analysis2/plots.MTA2.rda")
load("data analysis2/plots.regression.rda")
plt.si.nn=re.all%>%filter(si.nn.t12.por>0)

library(ggpubr);library(dplyr)
dat=MPA %>% select(PLT_CN,mpa.dead,mpa.native.survior) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mpa."),names_to = "Type", values_to = "MPA")#%>%filter(PLT_CN%in%plt.si.nn$PLT_CN)
t.test(MPA ~ Type,data= dat,paired = T)
dat2=MPA %>% select(PLT_CN,mpa.dead,mpa.native.survior) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mpa."),names_to = "Type", values_to = "MPA")%>%filter(PLT_CN%in%plt.si.nn$PLT_CN)
t.test(MPA ~ Type,data= dat2,paired = T)

plt1=dat%>%group_by(Type)%>%summarise(value=mean(MPA),se=sd(MPA)/sqrt(length(MPA)))
plt2=dat2%>%group_by(Type)%>%summarise(value=mean(MPA),se=sd(MPA)/sqrt(length(MPA)))
plt1$Type=factor(plt1$Type,levels=c("mpa.native.survior","mpa.dead"))
plt2$Type=factor(plt2$Type,levels=c("mpa.native.survior","mpa.dead"))

p1=ggplot(plt1) +
   geom_bar( aes(x=Type, y=value), stat="identity", fill="gray", alpha=0.7)+
	coord_cartesian(ylim=c(280,360))+
    geom_errorbar(aes(x=Type, ymin=value-1.96*se, ymax=value+1.96*se), width=0.4, colour="black", alpha=0.9, size=1.3)+
	labs(x="Type of species",y="Mean phylogenetic distance")+
	#scale_x_discrete(labels=c("mta.dead" = "Extinct natives vs.\nNative surviors","mta.native.survior" = "Native surviors vs.\nThemselves", "mta.si.nn" = "Non-native colonizers vs.\nNative surviors","mta.si.nn.dead"="Non-native colonizers vs.\nExtinct natives"))+
	scale_x_discrete(labels=c("mpa.dead" = "Extinct natives vs.\nNative surviors","mpa.native.survior" = "Native surviors vs.\nThemselves"))+	
		theme_classic()+theme(axis.title.x=element_blank(),
		axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=12,color='black'))+
	annotate("text", x = 1, y = 358, label = "a", size = 8, color = "darkgray")+
	annotate("text", x = 2, y = 358, label = "a", size = 8, color = "darkgray")
p2=ggplot(plt2) +
   geom_bar( aes(x=Type, y=value), stat="identity", fill="gray", alpha=0.7)+
	coord_cartesian(ylim=c(280,360))+
    geom_errorbar(aes(x=Type, ymin=value-1.96*se, ymax=value+1.96*se), width=0.4, colour="black", alpha=0.9, size=1.3)+
	labs(x="Type of species",y="Mean phylogenetic distance")+
	#scale_x_discrete(labels=c("mta.dead" = "Extinct natives vs.\nNative surviors","mta.native.survior" = "Native surviors vs.\nThemselves", "mta.si.nn" = "Non-native colonizers vs.\nNative surviors","mta.si.nn.dead"="Non-native colonizers vs.\nExtinct natives"))+
	scale_x_discrete(labels=c("mpa.dead" = "Extinct natives vs.\nNative surviors","mpa.native.survior" = "Native surviors vs.\nThemselves"))+	
		theme_classic()+theme(axis.title.x=element_blank(),
		axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=12,color='black'),
		axis.text.y = element_blank(),
		axis.title.y = element_blank())+
	annotate("text", x = 1, y = 358, label = "a", size = 8, color = "darkgray")+
	annotate("text", x = 2, y = 358, label = "a", size = 8, color = "darkgray")
ggarrange(p1,p2,ncol=2,align="h",labels="auto",hjust=-1,vjust=1.5)
