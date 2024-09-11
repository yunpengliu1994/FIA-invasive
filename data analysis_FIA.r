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
	state <-  readFIA(paste('data analysis2/FIA20240327/',statename,sep=""),tables=c('PLOT', 'TREE','SURVEY','COND','TREE_GRM_COMPONENT'))
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
			RECONCILECD,DIA,MICR_COMPONENT_AL_FOREST,SUBP_COMPONENT_AL_FOREST,CONDPROP_UNADJ,STATUSCD,COND_STATUS_CD)%>%
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
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
inv.plots=fia%>%filter(degreeOfEstablishment%in%status[1:2])%>%select(pltID)%>%distinct
inv.plots$invaded=1
fia=fia%>%left_join(inv.plots,by="pltID")
fia[is.na(fia$invaded),"invaded"]=0
fia[fia$Accepted_name=="Salix x pendulina nothof. x salamonii","Accepted_name"]="Salix x sepulcralis"
fia$DIA2=ifelse(fia$DIA>=12.5,"tree","sapling")
#fia$INVR = substr(fia$pltID,nchar(fia$pltID)-3,nchar(fia$pltID))

save(fia,file="data analysis2/fia.abundance.rda")

allplots=fia%>%filter(LON>=-100&MEASYEAR>=2000&MEASYEAR<=2021)
plt.pre=allplots%>%select(PLT_CN,PREV_PLT_CN)%>%distinct%>%na.omit
allplots.pre=fia%>%filter(PLT_CN%in%plt.pre$PREV_PLT_CN)
allplots.after=fia%>%filter(PREV_PLT_CN%in%plt.pre$PLT_CN)

invaded.plots=fia%>%filter(LON>=-100&invaded==1)
plt.pre=invaded.plots%>%select(PLT_CN,PREV_PLT_CN)%>%distinct%>%na.omit
invaded.pre=fia%>%filter(PLT_CN%in%plt.pre$PREV_PLT_CN)
invaded.after=fia%>%filter(PREV_PLT_CN%in%plt.pre$PLT_CN)

allplots=rbind(allplots.pre,allplots,allplots.after,invaded.plots,invaded.pre,invaded.after)%>%
	distinct%>%filter(!is.na(Abund_weight))
save(allplots,file="data analysis2/fia.abundance.ena.rda")

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
load("data analysis2/fia.abundance.ena.rda")
allplots=allplots%>%filter(MEASYEAR>=2000&MEASYEAR<=2021)
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
			legend.title=element_text(face="bold",size=15),
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
		legend.title=element_text(face="bold",size=15))	
		
######################################################
##richness and biomass changes through time #######
######################################################
library(dplyr);library(ggpubr);require(ggpmisc)
#detach(package:raster)
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
load("data analysis2/fia.abundance.ena.rda")
allplots=allplots%>%filter(LON>=-100&MEASYEAR>=2000&MEASYEAR<=2021)

get.stat=function(allplots.t,vars,ecoregion=FALSE){
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
			inv.plots.abd=allplots.t%>%filter(invaded==1)%>%group_by(MEASYEAR,pltID)%>%
					summarize(nat.abd=length(unique(VARS[which(is.na(degreeOfEstablishment))])),inv.abd=length(unique(VARS[which(degreeOfEstablishment%in%status[1:2])])),.groups="keep")
			uninv.plots.abd=allplots.t%>%filter(invaded==0&is.na(degreeOfEstablishment))%>%distinct%>%group_by(MEASYEAR,pltID)%>%summarize(n=length(unique(VARS)),.groups="keep")
		}else{
			inv.plots.abd=allplots.t%>%filter(invaded==1)%>%group_by(MEASYEAR,pltID)%>%
					summarize(nat.abd=sum(VARS[which(is.na(degreeOfEstablishment))]),inv.abd=sum(VARS[which(degreeOfEstablishment%in%status[1:2])]),.groups="keep")
			uninv.plots.abd=allplots.t%>%filter(invaded==0&is.na(degreeOfEstablishment))%>%distinct%>%group_by(MEASYEAR,pltID)%>%summarize(n=sum(VARS),.groups="keep")	
		}
		abd.inv=inv.plots.abd%>%group_by(MEASYEAR)%>%summarize(nPlots=length(pltID),native.mean=mean(nat.abd),inv.mean=mean(inv.abd),native.se=sd(nat.abd)/sqrt(length(nat.abd)),inv.se=sd(inv.abd)/sqrt(length(inv.abd)))
		abd.nat=uninv.plots.abd%>%group_by(MEASYEAR)%>%summarize(nPlots=length(pltID),univ.native.mean=mean(n),univ.native.se=sd(n)/sqrt(length(n)))
		native=abd.inv%>%select(MEASYEAR,nPlots,native.mean,native.se);inv=abd.inv%>%select(MEASYEAR,nPlots,inv.mean,inv.se)
		colnames(abd.nat)=colnames(native)=colnames(inv)=c("MEASYEAR","nPlots","abd.mean","abd.se")	
	}
	#mean richness of plots through time
	abd.stat=rbind(data.frame(Type0="Native",Type="Native(Invaded plots)",native),data.frame(Type0="Non-native",Type="Non-native",inv),data.frame(Type0="Native",Type="Native(Uninvaded plots)",abd.nat))
	abd.stat[is.na(abd.stat)]=0
	abd.stat$Type0=factor(abd.stat$Type0,levels=c("Non-native","Native"))
	abd.stat$Type=factor(abd.stat$Type,levels=c("Non-native","Native(Invaded plots)","Native(Uninvaded plots)"))
	return(abd.stat)
}
#ba.stat=get.stat(allplots,"Abund_weight")# basal area through time
bio.stat=get.stat(allplots,"Biomass")
sp.stat=get.stat(allplots,"Accepted_name")
theme=theme_bw()+
	theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),		
		strip.text = element_text(colour = 'black', size = rel(1.5)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))
p.sp=ggplot(sp.stat, aes(x=MEASYEAR, y=abd.mean,fill=Type,shape=Type)) + 	
		geom_errorbar(aes(ymin=abd.mean-1.96*abd.se, ymax=abd.mean+1.96*abd.se,color=Type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
		scale_fill_manual(values=c("#F8766D","#00BA38", "#619CFF"))+	
 		scale_shape_manual(values=c(21,21,24))+ 	
		geom_point(size=3,color="black",show.legend=T) + 		
		labs(x="Year",y="Mean species richness within plot") +
		geom_smooth(aes(x = MEASYEAR, y = abd.mean,color=Type),data=sp.stat,method="lm",linewidth=1.5,show.legend=FALSE,se =T,linetype=1)+		
		#stat_poly_eq(aes(label =  paste(stat(p.value.label), "*\"\"",sep = ""),color=Type),rr.digits = 2,p.digits = 2,label.x=0.1,label.y=c(1,0.92),size=6)+
		facet_wrap(~Type0, scales="free_y",ncol = 1)+
		theme+theme(legend.position= c(0.65,0.6),
			legend.background = element_rect(fill = NA),
			legend.title=element_blank(),
			legend.text=element_text(face="bold",size=15))	
p.bio=ggplot(bio.stat, aes(x=MEASYEAR, y=abd.mean,fill=Type,shape=Type)) + 	
		geom_errorbar(aes(ymin=abd.mean-1.96*abd.se, ymax=abd.mean+1.96*abd.se,color=Type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
		scale_fill_manual(values=c("#F8766D","#00BA38", "#619CFF"))+	
 		scale_shape_manual(values=c(21,21,24))+ 
		geom_point(size=3,color="black",show.legend=F) + 		
		labs(x="Year",y="Mean biomass within plot (Mg/Ha)") +
		geom_smooth(aes(x = MEASYEAR, y = abd.mean,color=Type),data=bio.stat,method="lm",linewidth=1.5,show.legend=FALSE,se =T,linetype=1)+
		#stat_poly_eq(aes(label =  paste(stat(p.value.label), "*\"\"",sep = ""),color=Type),rr.digits = 2,p.digits = 2,label.x=0.1,label.y=c(1,0.92),size=6)+			
		facet_wrap(~Type0, scales="free_y",ncol = 1)+
		theme
ggarrange(p.sp,p.bio,nrow=1,labels="auto")		

mf=function(div,dat.plt){
	dat=dat.plt%>%filter(Type==div)
	m=lm(abd.mean~MEASYEAR,dat)
	re=c(Type=div,df= summary(m)$df[2],F = round(summary(m)$f[1],2))
	return(re)
}
rbind(do.call(rbind,lapply(c("Native(Invaded plots)","Non-native","Native(Uninvaded plots)"),mf,sp.stat)), do.call(rbind,lapply(c("Native(Invaded plots)","Non-native","Native(Uninvaded plots)"),mf,bio.stat)))

## time trend in richness for each ecoregion ----
sp.stat=get.stat(allplots,"Accepted_name",ecoregion=TRUE)
#remove the ecoregion of a certain year when the number of invaded/uninvaded plots <=5 and then remove ecoregion do not contain invaded plots
sp.stat=sp.stat%>%filter(nPlots>5)
low.eco=sp.stat%>%select(Type0,ECOSUBCD)%>%distinct%>%group_by(ECOSUBCD)%>%tally%>%filter(n==2)
sp.stat=sp.stat%>%filter(ECOSUBCD%in%low.eco$ECOSUBCD)
#plot for each ecoregion
p=list();ecoregion=unique(sp.stat$ECOSUBCD)
for (i in ecoregion){
	dat=sp.stat%>%filter(ECOSUBCD==i)
	if (length(unique(dat$Type0))<2) next
	p[[i]]=ggplot(dat, aes(x=MEASYEAR, y=abd.mean,fill=Type,shape=Type)) + 	
		geom_errorbar(aes(ymin=abd.mean-1.96*abd.se, ymax=abd.mean+1.96*abd.se,color=Type),width=0.5,size=0.8,show.legend=FALSE,alpha=0.8)+
		scale_fill_manual(values=c("#F8766D","#00BA38", "#619CFF"))+	
 		scale_shape_manual(values=c(21,21,24))+ 	
		geom_point(size=2,color="black",show.legend=T) + 		
		labs(title=paste("Ecoregion",i),x="Year",y="Mean species richness within plot") +
		geom_smooth(aes(x = MEASYEAR, y = abd.mean,color=Type),data=dat,method="gam",linewidth=1.5,show.legend=FALSE,se =T,linetype=1)+		
		#stat_poly_eq(aes(label =  paste(stat(p.value.label), "*\"\"",sep = ""),color=Type),rr.digits = 2,p.digits = 2,label.x=0.1,label.y=c(1,0.92),size=3)+
		facet_wrap(~Type0, scales="free_y",ncol = 1)+
		theme_bw()+
		theme(axis.text = element_text(size=12,color='black'),
		axis.text.x = element_text(angle=15),
		axis.title = element_blank(),		
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(1), linetype = 2),
		legend.position= c(0.65,0.6),
		legend.background = element_rect(fill = NA),
		legend.title=element_blank(),
		legend.text=element_text(face="bold",size=12))+ 
		guides(fill = guide_legend(override.aes = list(size=3)))	
	if (which(ecoregion==i)<8) p[[i]]=p[[i]]+theme(axis.text.x = element_blank())
}
p1=ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],align="h",common.legend=T,legend="right")

mf=function(div,dat.plt){
	dat=dat.plt%>%filter(Type==div)
	m=lm(abd.mean~MEASYEAR,dat)
	re=data.frame(Type=div,p=round(summary(m)$coef[2,4],2) ,df= summary(m)$df[2],F = round(summary(m)$f[1],2))
	return(re)
}
sta.all=c();eco=unique(sp.stat$ECOSUBCD)
for (i in eco){
dat.eco=sp.stat%>%filter(ECOSUBCD%in%i)
if(nrow(dat.eco)>0) {sta=do.call(rbind,lapply(c("Native(Invaded plots)","Non-native","Native(Uninvaded plots)"),mf,dat.eco));sta$Ecoregion=i}
sta.all=rbind(sta.all,sta)
}
sta.all=sta.all%>% tidyr::pivot_wider(names_from = Type, values_from = c(p,df,F))
write.csv(sta.all,"data analysis2/statisticsLM.csv")

library(sf)
ecomap=st_read("S_USA.EcoMapProvinces/S_USA.EcoMapProvinces.shp")%>% st_transform("+proj=longlat +datum=WGS84")
nc3_points <- sf::st_point_on_surface(ecomap)
nc3_coords <- as.data.frame(sf::st_coordinates(nc3_points))
nc3_coords$MAP_UNIT_S <- ecomap$MAP_UNIT_S
nc3_coords[nc3_coords$MAP_UNIT_S=="Water","MAP_UNIT_S"]=""
#nc3_coords=nc3_coords%>%filter(MAP_UNIT_S%in%ecoregion)
p2=ggplot() + geom_sf(data=ecomap,color="darkgray",fill="black",show.legend=F) +				
		geom_text(data = nc3_coords, aes(X, Y, label = MAP_UNIT_S), colour = "white",size=4.5)+
		xlim(-102,-67)+theme_bw()+
		theme(panel.background = element_rect(fill = "#619CFF", colour = 'black'),
			axis.text = element_blank(),
			axis.title = element_blank(),
			panel.grid =element_blank())	
p1+patchwork::inset_element(p2, left = 0.61, bottom = 0.01, right = 0.84, top = 0.3)
ggarrange(p1,p2)		
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
load("data analysis2/fia.abundance.ena.rda")
load("data analysis2/traits.impute.rda")
load("data analysis2/phylo.imputate.Rdata")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
#remove established species
fia=allplots%>%filter(LON>=-100&(is.na(degreeOfEstablishment)|degreeOfEstablishment%in%status[1:2]))
fia$invasive=ifelse(is.na(fia$degreeOfEstablishment),"native","invasive")

#make abundance matrix
abund_mat <- fia %>%  group_by(pltID,Accepted_name)  %>% summarize(abund_w = sum(Abund_weight),.groups ="keep")%>% 
		tidyr::pivot_wider(names_from = Accepted_name, values_from = abund_w) %>% tibble::column_to_rownames("pltID")	
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
	colnames(FRic)[1]="pltID"	
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
		
## structural diversity ---
# ref: Relationships between stand growth and structural diversity in spruce-dominated forests in New Brunswick, Canada
# With the Shannon–Wiener index approach, DBH and height had
# to be grouped into discrete classes. For DBH, 2, 4, and
# 6 cm classes were tested, and for height, 2, 3, 4, and 5 m
# classes were tested to calculate the index.
fia$dia.class=round(fia$DIA/2)*2 #dia.class cm
#tapply(fia$pltID,fia$dia.class,length)
fia$dia.class=ifelse(fia$dia.class>74,">74",fia$dia.class)
fia$HT.class= round(fia$HT) #ht.class m
#tapply(fia$pltID,fia$HT.class,length)
fia$HT.class=ifelse(fia$HT.class>39,">39",fia$HT.class)
# Tree species diversity index: Shannon–Wiener index for species
Hs=fia %>% select(pltID,SPCD,Abund_weight) %>%group_by(pltID)%>% mutate(p = Abund_weight / sum(Abund_weight,na.rm =T))%>%
	group_by(pltID,SPCD) %>% summarize(p.i = sum(p,na.rm =T)) %>%  group_by(pltID) %>%
	summarize(Hs = -sum(p.i*log10(p.i)))
# Tree size diversity index: Shannon–Wiener index by diameter classes
Hd=fia %>% select(pltID,dia.class,Abund_weight) %>%group_by(pltID)%>% mutate(p = Abund_weight / sum(Abund_weight,na.rm =T))%>%
	group_by(pltID,dia.class) %>% summarize(p.i = sum(p,na.rm =T)) %>%  group_by(pltID) %>%
	summarize(Hd = -sum(p.i*log10(p.i)))
# Tree height diversity index: Shannon–Wiener index by height classes
Hh=fia %>% filter(!is.na(HT.class))%>% select(pltID,HT.class,Abund_weight) %>% 
	group_by(pltID) %>% mutate(p = Abund_weight / sum(Abund_weight,na.rm =T))%>%
	group_by(pltID,HT.class) %>% summarize(p.i = sum(p,na.rm =T)) %>%  group_by(pltID) %>%
	summarize(Hh = -sum(p.i*log10(p.i)))
# Mean structural diversity index: Mean value of tree species, size, and height indices	
Hsdh = Hs  %>%left_join(Hd,by="pltID") %>% left_join(Hh,by="pltID")
Hsdh$Hsdh=(Hsdh$Hs + Hsdh$Hd + Hsdh$Hh)/3
# Integrated diversity index of tree species and size: Integrated Shannon–Wiener index for species and diameter
Hsd=fia %>% select(pltID,SPCD,dia.class,Abund_weight) %>%group_by(pltID)%>% mutate(p = Abund_weight / sum(Abund_weight,na.rm =T))%>%
	group_by(pltID,SPCD,dia.class) %>% summarize(p.i = sum(p,na.rm =T)) %>%  group_by(pltID) %>%
	summarize(Hsd = -sum(p.i*log10(p.i)))
# Species profile index: Shannon–Wiener index calculation for the proportion of tree species in different
# stand layers; indicates integrated diversity of species and height	
Hsp = fia %>% filter(!is.na(HT)) %>%group_by(pltID)%>%
	summarize(Hmax.class = ifelse(HT>quantile(HT,0.8),1,ifelse(HT<=quantile(HT,0.5),3,2)),
		SPCD=SPCD,Abund_weight=Abund_weight) %>%
	group_by(pltID)%>% mutate(p = Abund_weight / sum(Abund_weight,na.rm =T))%>%
	group_by(pltID,SPCD,Hmax.class) %>% summarize(p.i = sum(p,na.rm =T)) %>%  group_by(pltID) %>%
	summarize(Hsp = -sum(p.i*log10(p.i)))
# Gini coefficient for DBH: Measurements of the deviation from perfect equality
# GCh <0 if DBHi>DBHj & HTi<HTj, e.g.,pltID==1_12_19_49_2006
# GCh = Na if only 1 tree are in the plot, e.g., pltID==1_17_77_20206_2010
GCd = fia %>% group_by(pltID) %>% 
	summarize(SPCD=SPCD,j = rank(DIA),n=length(DIA),Abund_weight=Abund_weight) %>% 
	summarize(numerator=(2*j-n-1)*Abund_weight,denominator=(n-1)*Abund_weight) %>% 
	group_by(pltID) %>% summarize(GCd = sum(numerator)/sum(denominator))
GCh = fia %>% filter(!is.na(HT)) %>% group_by(pltID) %>% 
	summarize(SPCD=SPCD,j=rank(HT,na.last="keep"),n=length(HT),Abund_weight=Abund_weight) %>% 
	summarize(numerator=(2*j-n-1)*Abund_weight,denominator=(n-1)*Abund_weight) %>% 
	group_by(pltID) %>% summarize(GCh = sum(numerator)/sum(denominator))

fin.inv=cbind(div1.inv,fin.inv);fin.nat=cbind(div1.nat,fin.nat)
fin=cbind(div1,fin)
colnames(fin.inv)[-7]=paste(colnames(fin.inv)[-7],"invasive",sep=".")
colnames(fin.nat)[-7]=paste(colnames(fin.nat)[-7],"native",sep=".")
div.all=list(fin,Hsdh,Hsd,Hsp,GCd,GCh,fin.inv,fin.nat)%>% purrr::reduce(left_join,by='pltID')

mf=function(x){
	stat=length(na.omit(x))
	if (stat>0){
		if (stat==1) return(x) 
		if (stat>1) {a=table(x);return(names(a)[which(a==max(a))][1])}
	}else {return(NA)}
}
plots.div=fia%>%group_by(pltID)%>%summarize(pltID.no.year=unique(pltID.no.year),year=unique(MEASYEAR),invaded=unique(invaded),
	ecoregion=mf(ECOSUBCD),soiltype=mf(PHYSCLCD),disturb=mean(DSTRBCD1,na.rm=T),
	biomass=mean(DRYBIO_AG+DRYBIO_BG,na.rm=T),std=mean(STDAGE,na.rm=T),BA=sum(Abund_weight))%>%
	left_join(div.all,by='pltID')
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
theme=theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),		
		strip.text = element_text(colour = 'black', size = rel(1.5)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2),
		legend.position= "right",
		legend.background = element_rect(fill = NA),
		legend.title=element_text(face="bold",size=15),
		legend.text=element_text(face="bold",size=12))
ggplot(dat.plt, aes(x=year, y=mean,fill=invaded)) + 		
	geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se,color=invaded),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+		
	geom_point(shape=21,size=2,color="black",alpha=0.8) + 		
	labs(x="Year",y="Mean diversity within plot",fill="Plot type") +
	geom_smooth(aes(x = year, y = mean,color=invaded),data=dat.plt,method="lm",linewidth=1,show.legend=FALSE,se =T,linetype=1,alpha=0.5)+
	#stat_poly_eq(aes(label =  paste(stat(p.value.label), "*\"\"",sep = ""),color=invaded),rr.digits = 2,p.digits = 2,label.x=1,label.y=c(1,0.92),size=5)+			
	scale_fill_manual(values=c("#F8766D","#619CFF"),labels=c("Invaded plots","Uninvaded plots"))+ 	
	scale_color_manual(values=c("#F8766D","#619CFF"),labels=c("Invaded plots","Uninvaded plots"))+ 		
	facet_wrap(~index, scales="free_y",ncol=2,labeller = labeller(index = index.labs))+			
	theme+ guides(fill = guide_legend(override.aes = list(size=5)))

mf=function(div,dat.plt,inv){
	dat=dat.plt%>%filter(index==div&invaded==inv)
	m=lm(mean~year,dat)
	re=c(div=div,inv=inv,summary(m)$coef[2,],Adj.r = summary(m)$adj, F = summary(m)$f[1])
	return(re)
}
rbind(do.call(rbind,lapply(c("FAD.native","M_Ap.native","M_At.native","PD.native"),mf,dat.plt,inv=1)), do.call(rbind,lapply(c("FAD.native","M_Ap.native","M_At.native","PD.native"),mf,dat.plt,inv=0)))


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
#plots=allplots%>%select(PLT_CN,pltID.no.year)%>%distinct%>%group_by(pltID.no.year)%>%tally%>%filter(n>2)

get.plt.reg=function(plt,allplots,sp.status,tree,trait_mat){
#re.all=c()
#for(plt in plots.cn$PLT_CN[which(plots.cn$PLT_CN==plt):nrow(plots.cn)]){
# plt=plots.cn$PLT_CN[1230]
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
save(re.all,file="/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.regression2.rda")

library(dplyr)
load("data analysis2/plots.regression.rda")
# hist(re.all$measure.time12,breaks=20);abline(v=4,col="red");abline(v=6,col="red");abline(v=3,col="blue");abline(v=8,col="blue");
# windows()
# hist(re.all$measure.time23,breaks=20);abline(v=4,col="red");abline(v=6,col="red");abline(v=3,col="blue");abline(v=8,col="blue");
# plot(dat$si.tot.t12,dat$inv.tot.t12);abline(a=0,b=1,col="red")
# dat%>%filter(si.tot.t12==0&inv.tot.t12>10)%>%nrow

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
		
		xvars.re=c("Bcol.non_native","Bsur.non_native","Bcol.native","Biomass.t2","Diversity.t2","STDAGE","measure.time12","measure.time23")
		xvar.labs <- c("Bcol.non_native","Bsur.non_native","Bcol.native","B2","Y2","StandAge","t12","t23")
		
		formu=as.formula(paste("yvar ~",paste(xvars.re,collapse="+"),"+(1|ECOSUBCD)",sep=" "))				
	}else{
		# xvars.re=c("si.nn.t12.por","nn.t2","suvor.t1.biomass","si.tot.t12","Biomass.t2","STDAGE","measure.time12","measure.time23")
		# xvar.labs <- c("Colonizer.Non-native", "Non-native","Survivor","Colonizer","Biomass","StandAge","t12","t23")
		xvars.re=c("Bcol.non_native","Bsur.non_native","Bcol.native","Biomass.t2","STDAGE","measure.time12","measure.time23")
		xvar.labs <- c("Bcol.non_native","Bsur.non_native","Bcol.native","B2","StandAge","t12","t23")
		
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
stat.all2=stat.all%>%select(-type)%>%tidyr::pivot_wider(names_from = var, values_from = c(Estimate,se,p))%>%
	select(yvar,contains(xvar.labs[1]),contains(xvar.labs[2]),contains(xvar.labs[3]),contains(xvar.labs[4]),contains(xvar.labs[5]),contains(xvar.labs[6]),contains(xvar.labs[7]),contains(xvar.labs[8]))	
# write.csv(stat.all2,"data analysis2/Fig.3.data.csv")
stat.all$yvar=factor(stat.all$yvar,levels=yvars)
ggplot() +geom_errorbar(data=stat.all,aes(x=Estimate, y=var,xmin=Estimate-1.96*se, xmax=Estimate+1.96*se,color=type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
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
	
######################################################
#are individuals that disappeared in t3 more likely to be small and rare?
######################################################
# Compare the dead trees in t3 among each case in their num. of individuals of the species, DIA, HT, BA, biomass and their differences with the biggese in the PLOT
library(dplyr)
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/fia.abundance.rda")
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.regression.rda")
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/traits.impute.rda")

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
		
	dorminant=t2%>%filter(STATUSCD==1)%>%group_by(Accepted_name)%>%summarise(Abd=length(Accepted_name),Biomass=sum(Biomass,na.rm = T))
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
no_cores <- 27#detectCores() - 1
mycl <- makePSOCKcluster(no_cores); 
re.all.dead=do.call(rbind,parLapply(cl=mycl,re.all$PLT_CN,get.dead.rare,allplots,trait_mat))
stopCluster(mycl)
save(re.all.dead,file="/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.dead.trees.V2.rda")

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
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/fia.abundance.rda")
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.regression.rda")
load("/blue/matthewthomas1/yunpeng.liu/data analysis2/traits.impute.rda")

# load("data analysis2/traits.impute.rda")
# load("data analysis2/fia.abundance.rda")
# load("data analysis2/plots.regression.rda")
trait_mat[is.na(trait_mat)]=0
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=fia%>%filter(STATUSCD%in%c(1,2)&LON>=-100&CONDPROP_UNADJ==1&(!degreeOfEstablishment%in%status[3]))#STATUSCD==1
allplots$measure.time=allplots$MEASYEAR+allplots$MEASMON/12
sp.status=allplots%>%filter(!is.na(PREV_PLT_CN))%>%select(Accepted_name,degreeOfEstablishment)%>%distinct
sp.status$inv=ifelse(is.na(sp.status$degreeOfEstablishment),0,1)

get.mta=function(plt,trait_mat,allplots,sp.status){
	#plt=re.all$PLT_CN[1]
	require(dplyr);require(FD)
	t2=allplots%>%filter(PLT_CN%in%plt)
	t3=allplots%>%filter(PREV_PLT_CN%in%unique(t2$PLT_CN))
	dead.t3=t3%>%filter(STATUSCD==2&(!is.na(PREV_TRE_CN))) #tree dead in t3, including 1) tree live in t2 while dead in t3; 2)tree dead in t2 and still there in t3
	if (nrow(dead.t3)==0) return(NULL)
	dead.t2=t2%>%filter(CN%in%dead.t3$PREV_TRE_CN&STATUSCD==1&is.na(degreeOfEstablishment))#native tree live in t2 while dead in t3
	if (nrow(dead.t2)==0) return(NULL)	
	new.inv.t2=t2%>%filter(is.na(PREV_TRE_CN))%>%distinct
	if (nrow(new.inv.t2)==0) return(NULL) #next; ## No new individuals in t2
	t1=allplots%>%filter(PLT_CN%in%unique(t2$PREV_PLT_CN))	
	if (nrow(t1)==0) return(NULL) #next; ## No plot mears or not used in t1
	sp.all=unique(t2$Accepted_name)
	if (length(sp.all)<2)	return(NULL)
	
	new.sp.suvor=new.inv.t2%>%select(Accepted_name,CN,Biomass)%>%left_join(sp.status,by="Accepted_name")
	new.sp.suvor$suvor.t1=new.sp.suvor$Accepted_name%in%t1$Accepted_name #t2 new individuals of which species was present at t1
	new.sp.suvor$suvor.t3=new.sp.suvor$CN%in%t3$PREV_TRE_CN #t2 new individuals of which species survived to t3 
		
	si.native=unique(new.sp.suvor[new.sp.suvor$suvor.t1==F&new.sp.suvor$suvor.t3==T&new.sp.suvor$inv==0,"Accepted_name"])#successful non-native invader and the species was not present at t1.
	si.nn=unique(new.sp.suvor[new.sp.suvor$suvor.t1==F&new.sp.suvor$suvor.t3==T&new.sp.suvor$inv==1,"Accepted_name"])#successful non-native invader and the species was not present at t1.
	nn=sp.all[sp.all%in%sp.status[sp.status$inv==1,"Accepted_name"]]
	dead.sp=unique(dead.t2$Accepted_name)
	native.survior=t2%>%filter(STATUSCD==1&(!CN%in%dead.t2)&(!Accepted_name%in%sp.status[sp.status$inv==1,"Accepted_name"]))%>%as.data.frame%>%.[,"Accepted_name"]%>%unique#native species survived to t3
	## caculate Abundance weighted functional distinctiveness
	abd=t2%>%group_by(Accepted_name)%>%tally
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
		if (length(j.sp)==0) {
			return(NA)
		}else{
			dij=as.matrix(G.trait.mat[j.sp,i.sp]);colnames(dij)=i.sp
			nj=abd$n[abd$Accepted_name%in%j.sp];names(nj)=j.sp
			MTA=sum(dij*nj)/sum(nj)
			return(MTA)
		}		
	}
	re=data.frame(PLT_CN=plt,measure.time12=unique(t2$measure.time)-unique(t1$measure.time),
		measure.time23=unique(t3$measure.time)-unique(t2$measure.time),
		mta.nn=ifelse(length(nn)==0,NA,do.call(mean,lapply(nn,mta.cal,G.trait.mat,abd,j.sp=NULL,sp.all=sp.all))),nn=length(nn),
		mta.si.native=ifelse(length(si.native)==0,NA,do.call(mean,lapply(si.native,mta.cal,G.trait.mat,abd,j.sp=NULL,sp.all=sp.all))),si.native=length(si.native),
		mta.si.nn=ifelse(length(si.nn)==0,NA,do.call(mean,lapply(si.nn,mta.cal,G.trait.mat,abd,j.sp=NULL,sp.all=sp.all))),si.nn=length(si.nn),	
		mta.dead=ifelse(length(dead.sp)==0,NA,do.call(mean,lapply(dead.sp,mta.cal,G.trait.mat,abd,j.sp=NULL,sp.all=sp.all))),dead=length(dead.sp),
		mta.native.survior=ifelse(length(native.survior)==0,NA,do.call(mean,lapply(native.survior,mta.cal,G.trait.mat,abd,j.sp=NULL,sp.all=sp.all))),native.survior=length(native.survior)	
	)	
	return(re)
}

library(parallel)
no_cores <- 37#detectCores() - 1
mycl <- makePSOCKcluster(no_cores); 
MTA=do.call(rbind,parLapply(cl=mycl,re.all$PLT_CN,get.mta,trait_mat,allplots,sp.status))
stopCluster(mycl)
save(MTA,file="/blue/matthewthomas1/yunpeng.liu/data analysis2/plots.MTA.rda")

# plot the Difference
# ref:http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
load("data analysis2/plots.MTA.rda")
load("data analysis2/plots.regression.rda")
library(ggpubr);library(dplyr)
dat1=MTA %>% select(PLT_CN,mta.native.survior,mta.si.nn) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")
dat2=MTA %>% select(PLT_CN,mta.native.survior,mta.si.native) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")
dat3=MTA %>% select(PLT_CN,mta.dead,mta.native.survior) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")
dat4=MTA %>% select(PLT_CN,mta.dead,mta.si.nn) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")
dat5=MTA %>% select(PLT_CN,mta.dead,mta.si.native) %>%na.omit%>% tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")

plt=re.all%>%filter(nn.t2!=0)
dat=dat3%>%filter(PLT_CN%in%plt$PLT_CN)
t.test(MTA ~ Type,data= dat,paired = T)
# re=compare_means(MTA ~ Type,data= dat, method = "t.test", paired = T)
# dat%>%group_by(Type)%>%summarize(MTA=mean(MTA),nplots=length(PLT_CN))

# dat=MTA %>% select(PLT_CN,mta.native.survior,mta.dead,mta.si.nn,mta.si.native) %>% na.omit()%>%tidyr::pivot_longer(cols=contains("mta."),names_to = "Type", values_to = "MTA")
# compare_means(MTA ~ Type,data= dat, method = "wilcox.test", paired = T)
# dat%>%group_by(Type)%>%summarize(MTA=mean(MTA),nplots=length(PLT_CN))

# ggboxplot(dat, x = "Type", y = "MTA")+stat_compare_means(paired = TRUE,label =  "p.signif", label.x = 1.5)

############################################################################################################################################################################################################################
##### code graveyard #######################################################################################################################################################################################################
############################################################################################################################################################################################################################

#conduct regression
library(lme4);
library(afex)#This will automatically add a p-value column to the output of the lmer for the fixed effects.
library(ggpubr)
#load("data analysis2/plots.dead.trees.rda")
load("data analysis2/plots.dead.trees.V2.rda")
load("data analysis2/plots.regression.rda")
#yvars=c("mta.dead","rare.individual.por.ba","rare.sp.por.ba","rare.Biomass.por.ba","rare.BA.por.ba","rare.individual.por","rare.sp.por","rare.Biomass.por","rare.BA.por")
#yvars=c("rare.individual.por","rare.sp.por","rare.Biomass.por")
yvars=c("mta.dead")
yvar.lab=c("FDis.loss")
names(yvar.lab)=yvars

re.all$Bcol.non_native=ifelse(is.na(re.all$si.nn.t12.por),0,re.all$si.nn.t12.por*re.all$si.tot.t12)
re.all$Bcol.native=ifelse(is.na(re.all$si.nn.t12.por),0,re.all$si.tot.t12-re.all$Bcol.non_native)
re.all$Bsur.non_native=ifelse(is.na(re.all$si.nn.t12.por),re.all$nn.t2,re.all$nn.t2-re.all$Bcol.non_native)

# xvars=c("si.nn.t12.por","nn.t2","suvor.t1.biomass","si.tot.t12","Mortality.rate","SR.t2","Biomass.t2","STDAGE","measure.time12","measure.time23")
# xvar.labs <- c("Colonizer.Non-native", "Non-native","Survivor","Colonizer","Mortality","SR","Biomass","StandAge","t12","t23")

xvars=c("SR.loss","Bcol.non_native","Bsur.non_native","Bcol.native","Biomass.t2","STDAGE","measure.time12","measure.time23")
xvar.labs <- c("SR.loss","Bcol.non_native","Bsur.non_native","Bcol.native","B2","StandAge","t12","t23")

formu=as.formula(paste("yvar ~",paste(xvars,collapse="+"),"+(1|ECOSUBCD)",sep=" "))
stat.all=c()
for(i in 1:length(yvars)){
	yvar=yvars[i];
	dat=re.all%>%left_join(re.all.dead,by="PLT_CN")%>%#filter(measure.time23>3&measure.time23<10&measure.time12>3&measure.time12<10)%>%#%>%filter(si.nn.t12.por==0&nn.t2==0)
		filter(!(Bcol.non_native==0&Bcol.native==0&Bsur.non_native==0))
	colnames(dat)[colnames(dat)==yvar]="yvar"
	dat[,xvars]=scale(dat[,xvars])
	M <- lmer(formu,data=dat)
	stat=as.data.frame(summary(M)$coef[-1,c(1,2,5)])
	colnames(stat)[2:3]=c("se","p")
	#rownames(stat)[rownames(stat)=="Diversity.t2"]=xvar
	stat$var=xvar.labs
	stat$var=factor(stat$var,levels=rev(xvar.labs))
	stat$type=ifelse((stat$Estimate+1.96*stat$se)<0,"negative",ifelse((stat$Estimate-1.96*stat$se)>0,"positive","insig"))
	stat$yvar=yvar#.lab
	stat.all=rbind(stat.all,stat)
}
stat.all$yvar=factor(stat.all$yvar,levels=yvars)
stat.all$var=factor(stat.all$var,levels=rev(xvar.labs))
p1=ggplot() +geom_errorbar(data=stat.all,aes(x=Estimate, y=var,xmin=Estimate-1.96*se, xmax=Estimate+1.96*se,color=type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
	geom_point(data=stat.all,aes(x=Estimate, y=var,fill=type),size=2,shape=21,color="black",show.legend=F) + 		
	scale_fill_manual(values=c("insig"="darkgray", "positive"="#6D9EC1","negative"="#E46726"))+	
	scale_color_manual(values=c("insig"="darkgray", "positive"="#6D9EC1","negative"="#E46726"))+	
	labs(x="Estimates") +
	# geom_text(data = stat.n%>%filter(inv=="native"), mapping = aes(x=0 , y=-0.5, label = paste0("Plots (Native) = ",n)), hjust   = "center", vjust   = "inward",size=4)+
	# geom_text(data = stat.n%>%filter(inv=="alien"), mapping = aes(x=0 , y=-1.2, label = paste0("Plots (NonNative) = ",n)), hjust   = "center", vjust   = "inward",size=4)+	
	facet_wrap(~yvar, scales="free_x",labeller = labeller(yvar=yvar.lab))+		
	theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),
		axis.title.y=element_blank(),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))+
	geom_vline(xintercept=0,col="lightgray",size=0.8,linetype="longdash")

dat=re.all%>%left_join(re.all.dead,by="PLT_CN")%>%#filter(measure.time23>3&measure.time23<10&measure.time12>3&measure.time12<10)%>%#%>%filter(si.nn.t12.por==0&nn.t2==0)
		#filter(!(Bcol.non_native==0&Bcol.native==0&Bsur.non_native==0))%>%
		select(SR.loss,mta.dead,rare.sp.por)%>%na.omit
colnames(dat)=yvar.lab
dat=dat%>%tidyr::pivot_longer(cols=SR.loss:PSR.loss.Rare,names_to = "Vars", values_to = "value")
dat$Vars=factor(dat$Vars,levels=yvar.lab)
p2=ggplot(dat, aes(x=value)) + 
  geom_histogram(color="black", fill="lightblue",bins=50)+
   labs(y="Number of Plots")+
  facet_wrap(~Vars, scales="free_x")+	
  theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),
		axis.title.x = element_blank(),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))
ggarrange(p1,p2,ncol=1,heights=c(1,0.5),align="v",labels="auto",hjust=0,vjust=0.8)

## triat space and phylogeny of non-natives to overall species 
#PCA 
load("data analysis2/traits.impute.rda")
load("data analysis2/fia.abundance.rda")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
splist=fia%>%filter(STATUSCD%in%c(1,2)&LON>=-100&CONDPROP_UNADJ==1&(!degreeOfEstablishment%in%status[3]))%>%select(Accepted_name,degreeOfEstablishment)%>%distinct
splist$invasive=ifelse(is.na(splist$degreeOfEstablishment),"Native","Non-native")
library(ggfortify)
trait_mat[is.na(trait_mat)]=0
trait_mat=data.frame(trait_mat,Accepted_name=rownames(trait_mat))%>%left_join(splist,by="Accepted_name")%>%filter(!is.na(invasive))
trait_mat$invasive=factor(trait_mat$invasive,levels=c("Non-native","Native"))
colnames(trait_mat)[1:5]=c("SLA","LN","Rmax","Hmax","WD")

pca_res <- prcomp(trait_mat[,c(1:5)], scale. = TRUE)
p1=autoplot(pca_res, data = trait_mat, fill = 'invasive',size=3,alpha=1,shape=21,color="black",
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 5, loadings.label.color ="black")+
scale_fill_manual(values=c("#F8766D","#619CFF"))+
geom_vline(xintercept=0,col="darkgray",size=1,linetype="longdash",alpha=0.5)+
geom_hline(yintercept=0,col="darkgray",size=1,linetype="longdash",alpha=0.5)+
guides(color=guide_legend(override.aes = list(size = 3)))+	
theme_bw()+theme(axis.text = element_text(size=15,color='black'),
	axis.title = element_text(face="bold",size=18,color='black'),
	legend.position=c(.82, .93),
	legend.background = element_rect(fill = NA,color="black"),
		legend.text=element_text(face="bold",size=15),
		legend.title=element_blank())

# if (!requireNamespace("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# BiocManager::install("ggtree")
library(ggtree);library(ape)
load("data analysis2/phylo.imputate.Rdata")
tre.plt=keep.tip(tree, tip=splist$Accepted_name[splist$Accepted_name%in%tree$tip])
splist$invasive=factor(splist$invasive,levels=c("Non-native","Native"))
p2 <- ggtree(tre.plt,layout="fan")  %<+%  splist + geom_tippoint(aes(color=invasive),size=2,show.legend=F)+scale_color_manual(values=c("#F8766D","transparent"))

ggarrange(p1,p2,labels="auto")

#all the possible cases of non-natives in each plot: Number of plot mears and the biomass and richness differences
# CASE 1:si.nn.t12.por>0&nn.t2>0 #(non-natives comes from both invaders and t1 survivors) : 228
# CASE 2:si.nn.t12.por==0&nn.t2>0 #(non-natives only comes from t1 survivors): 371
# CASE 3:si.nn.t12.por==0&nn.t2==0 #(no non-natives invaded in t2 or survived from t1, maybe this is the case for native invaders?): 14645
# CASE 4:is.na(si.nn.t12.por)&nn.t2>0 #(non-natives only comes from t1 survivors and all invaders in t2 are species present at t1): 533
# CASE 5:is.na(si.nn.t12.por)&nn.t2==0 #(no non-natives invaded in t2 or survived from t1, and all native invaders in t2 are species present at t1):27572

## statics on each cases:
dat=re.all%>%filter(measure.time23>3&measure.time23<10&measure.time12>3&measure.time12<10)#From 47152 to 26871
yvars=c("SR.t23","PD.t23","FAD.t23","M_Ap.t23","M_At.t23","Biomass.t23")
yvar.labs <- c("ΔSR", "ΔPD","ΔFAD","ΔMBL","ΔMTD","ΔBiomass")
names(yvar.labs) <- yvars
xvars=c("SR.t2","PD","FAD","M_Ap","M_Ap",NA)
dat2=dat%>%filter(si.nn.t12.por == 0 & nn.t2 == 0)
stat.all=c()
for(i in 1:length(yvars)){
	yvar=yvars[i];xvar=xvars[i]
	dat.t=dat2
	colnames(dat.t)[colnames(dat.t)==yvar]="yvar"
	if (!is.na(xvar)) {
		colnames(dat.t)[colnames(dat.t)==xvar]="Diversity.t2"
		xvars.re=c("si.tot.t12","suvor.t1.biomass","Diversity.t2","Biomass.t2","STDAGE","measure.time12","measure.time23")
		xvar.labs <- c("Colonizer","Survivor","Diveristy","Biomass","StandAge","t12","t23")		
		formu=as.formula(paste("yvar ~",paste(xvars.re,collapse="+"),"+(1|ECOSUBCD)",sep=" "))				
	}else{
		xvars.re=c("si.tot.t12","suvor.t1.biomass","Biomass.t2","STDAGE","measure.time12","measure.time23")		
		xvar.labs <- c("Colonizer","Survivor","Biomass","StandAge","t12","t23")
		formu=as.formula(paste("yvar ~",paste(xvars.re,collapse="+"),"+(1|ECOSUBCD)",sep=" "))		
	}
	dat.t[,xvars.re]=scale(dat.t[,xvars.re])
	M <- lmer(formu,data=dat.t)
	stat=as.data.frame(summary(M)$coef[-1,c(1,2,5)])
	colnames(stat)[2:3]=c("se","p")
	#rownames(stat)[rownames(stat)=="Diversity.t2"]=xvar
	stat$var=xvar.labs
	stat$var=factor(stat$var,levels=rev(xvar.labs))
	stat$type=ifelse((stat$Estimate+1.96*stat$se)<0,"negative",ifelse((stat$Estimate-1.96*stat$se)>0,"positive","insig"))
	stat$yvar=yvar
	stat.all=rbind(stat.all,stat)	
}

stat.all$yvar=factor(stat.all$yvar,levels=yvars)
ggplot() +geom_errorbar(data=stat.all,aes(x=Estimate, y=var,xmin=Estimate-1.96*se, xmax=Estimate+1.96*se,color=type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
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

##################################################
##### plots invaded between time 1 and 2 #########
##################################################
###   section 1 Changes of diversity and biomass in plots after being invaded  ----------
library(dplyr);library(ggpubr)
load("data analysis2/fia.abundance.rda")
load("data analysis2/plots.summary.rda")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
inv.sp=fia%>%filter(degreeOfEstablishment%in%status[1:2])%>%select(Accepted_name)%>%distinct%>%.[,"Accepted_name"]
inv.plots=fia%>%filter(Accepted_name%in%inv.sp)%>%select(pltID)%>%distinct
inv.plots$invaded=1
allplots=fia%>%filter((!is.na(Abund_weight))&LON>=-100&(is.na(degreeOfEstablishment)|degreeOfEstablishment%in%status[1:2]))
allplots=allplots%>%select(-invaded)%>%left_join(inv.plots,by="pltID")
allplots[is.na(allplots$invaded),"invaded"]=0
plots.div=plots.div%>%select(-invaded)%>%left_join(inv.plots,by="pltID")
plots.div[is.na(plots.div$invaded),"invaded"]=0
vars=c("HillEven.native","PD.native","M_Ap.native","FAD.native","M_At.native")
	
get.re=function(plots.div.null,null.sp,allplots.null,vars,rep=0){
	plt.inv=plots.div.null%>%filter(invaded==1)%>%select(pltID.no.year)%>%distinct #2430 plots
	plt.inv.multi=plots.div.null%>%select(year,pltID.no.year)%>%distinct%>%group_by(pltID.no.year)%>%tally%>%filter(n>1&pltID.no.year%in%plt.inv$pltID.no.year)#1931 plots ranges from 2 to 5 measurements
	plt.inv.alltime=plots.div.null%>%filter(pltID.no.year%in%plt.inv.multi$pltID.no.year)%>%select(year,pltID.no.year,invaded)%>%distinct%>%group_by(pltID.no.year)%>%
		summarize(invade.type=length(invaded)-sum(invaded))%>%#invade.type:0 invaded in all measurements, 805 out of 1944 plots; <0 invaded in some of the measurements, 1127 out of 1932 plots
		filter(invade.type>0)
	#vars=c("SR","BA","biomass","HillEven","shannon","PD","M_Ap","phy.shannon","FAD","M_At","Hsp")	
	allplots.inv=plots.div.null%>%filter(pltID.no.year%in%plt.inv.alltime$pltID.no.year)%>%
		select(pltID:disturb,invaded,all_of(vars))%>%
		group_by(pltID.no.year,invaded)%>%summarise(across(HillEven.native:M_At.native, ~ mean(.x, na.rm = TRUE)),.groups="keep")
	SR.native=allplots.null%>%filter(pltID.no.year%in%plt.inv.alltime$pltID.no.year&(!Accepted_name%in%null.sp))%>%group_by(pltID.no.year,invaded)%>%summarise(SR.native=length(unique(Accepted_name)),.groups="keep")
	BA.native=allplots.null%>%filter(pltID.no.year%in%plt.inv.alltime$pltID.no.year&(!Accepted_name%in%null.sp))%>%group_by(pltID.no.year,invaded)%>%summarise(BA.native=sum(Abund_weight),.groups="keep")
	biomass.native=allplots.null%>%filter(pltID.no.year%in%plt.inv.alltime$pltID.no.year&(!Accepted_name%in%null.sp))%>%group_by(pltID.no.year,invaded)%>%summarise(Biomass.native=sum(Biomass),.groups="keep")
	allplots.inv=list(SR.native,BA.native,biomass.native,allplots.inv)%>% purrr::reduce(left_join,by=c("pltID.no.year","invaded"))
	#add invaded plots where native species extinct after Invaded
	a=allplots.inv%>%filter(invaded==1);b=allplots.inv%>%filter(invaded==0)
	un=b$pltID.no.year[!b$pltID.no.year%in%a$pltID.no.year]
	if (length(un)>0){
		add=allplots.inv%>%filter(pltID.no.year%in%un)%>%as.data.frame
		add=replace(add,3:ncol(add),0);add=replace(add,2,1)
		allplots.inv=rbind(allplots.inv,add)
	}

	re=c()
	for (var.t in c("SR.native","BA.native","Biomass.native",vars)){
		dat=na.omit(allplots.inv[,c(var.t,"invaded","pltID.no.year")])
		dat[,1]=scale(dat[,1]);colnames(dat)[1]="var.t"	
		a=dat%>%filter(invaded==1);b=dat%>%filter(invaded==0)
		keep=intersect(b$pltID.no.year,a$pltID.no.year)#native species distinct after Invaded
		dat=dat%>%filter(pltID.no.year%in%keep)
		dat$invaded=factor(dat$invaded,levels=c(1,0))
		
		test=t.test(var.t ~ invaded,data=dat,alternative ="two.sided",paired = TRUE, conf.level = 0.95)
		P=test$p.value
		Px=ifelse(P<=0.001,"***",
			ifelse(P>0.001&P<=0.01,"**",
			ifelse(P>0.01&P<=0.05,"*",
			#ifelse(P>0.05&P<0.1,"",
			ifelse(P>0.05,"ns",P))))
		tmp=data.frame(vars=var.t,mean=test$est,se=test$std, P=P,p=ifelse(P<0.001,paste0("<0.001",Px),paste0(round(P,3),Px)))
		re=rbind(re,tmp)
	}
	re$type=ifelse(re$P>0.05,"insig",ifelse(re$mean>0,"positive","negative"))
	re$vars=factor(re$vars,levels=rev(c("SR.native","BA.native","Biomass.native",vars)))
	re$nplots=nrow(plt.inv.alltime)
	if(rep>0) re$rep=i
	return(re)
}

#plots invaded between time 1 and 2, and changes in invasive/native richness/abundance during time 2 and 3
re.inv=get.re(plots.div,inv.sp,allplots,vars)
p1=ggplot(re.inv) + 	
	geom_errorbar(aes(x=mean, y=vars,xmin=mean-1.96*se, xmax=mean+1.96*se,color=type),width=0.2,size=1,show.legend=FALSE,alpha=0.8)+
	geom_point(aes(x=mean, y=vars,fill=type),size=4,shape=21,color="black",show.legend=F) + 	
	scale_fill_manual(values=c("insig"="gray", "positive"="#6D9EC1","negative"="#E46726"))+	
	scale_color_manual(values=c("insig"="gray", "positive"="#6D9EC1","negative"="#E46726"))+		
	labs(x="Difference after invaded (scaled)") +	
	annotate("text",x=0.05 , y= 1,size=6,label = paste0("nPlots = ",unique(re.inv$nplots)),color = "black")+
	geom_vline(xintercept=0,col="darkgray",size=1,linetype="longdash")+
	theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),
		axis.title.y=element_blank())
#map these plots
library(sf)
ecoregion=st_read("S_USA.EcoMapProvinces/S_USA.EcoMapProvinces.shp")%>% st_transform("+proj=longlat +datum=WGS84")
plots.sp=allplots%>%filter(pltID.no.year%in%plt.inv.alltime$pltID.no.year)%>% dplyr::select(pltID.no.year,LON,LAT)%>% distinct%>%
	st_as_sf(coords=c("LON","LAT"))%>% st_set_crs("+proj=longlat +datum=WGS84")
p2=ggplot() + geom_sf(data=ecoregion,color="lightgray",fill ="white") +xlim(-100,-67)+ 
		geom_sf(data=plots.sp,color="#F8766D",size=1,show.legend=F)+
		theme_bw()+
		theme(axis.text = element_text(size=12,color='black'),			
			axis.title = element_blank(),
			panel.background = element_rect(fill = '#619CFF', colour = 'black'))	
ggarrange(p1,p2,ncol=2,align = "h",labels="auto")

## random sample the same amount of native species and then conduct the same comparison
native.sp=allplots%>%filter(is.na(degreeOfEstablishment))%>%select(Accepted_name)%>%distinct
set.seed(123)
re.all=c()
for (i in 1:30){
	null.sp=sample(native.sp$Accepted_name, length(inv.sp), replace = FALSE)
	null.plots=fia%>%filter(Accepted_name%in%null.sp)%>%select(pltID)%>%distinct
	null.plots$invaded=1
	allplots.null=allplots%>%filter(invaded==0)%>%select(-invaded)%>%left_join(null.plots,by="pltID")
	allplots.null[is.na(allplots.null$invaded),"invaded"]=0
	plots.div.null=plots.div%>%filter(invaded==0)%>%select(-invaded)%>%left_join(null.plots,by="pltID")
	plots.div.null[is.na(plots.div.null$invaded),"invaded"]=0
	re=get.re(plots.div.null,null.sp,allplots.null,vars,rep=i)	
	re.all=rbind(re.all,re)
	print(i)
}
save(re.all,file="data analysis2/null.effect.rda")
nplots=re.all%>%select(rep,nplots)%>%distinct%>%as.data.frame
ggplot(re.all) + 	
	geom_errorbar(aes(x=mean, y=vars,xmin=mean-1.96*se, xmax=mean+1.96*se,color=type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
	geom_point(aes(x=mean, y=vars,fill=type),size=2,shape=21,color="black",show.legend=F) + 	
	scale_fill_manual(values=c("insig"="gray", "positive"="#6D9EC1","negative"="#E46726"))+	
	scale_color_manual(values=c("insig"="gray", "positive"="#6D9EC1","negative"="#E46726"))+		
	labs(x="Difference after invaded (scaled)") +	
	#annotate("text",x=0.05 , y= 1,size=6,label = paste0("Plots = ",nrow(plt.inv.alltime)),color = "black")+
	geom_text(data = nplots, mapping = aes(x=0 , y=-0.5, label = paste0("nPlots = ",nplots)), hjust   = "center", vjust   = "inward",size=4)+	
	geom_vline(xintercept=0,col="darkgray",size=1,linetype="longdash")+
	facet_wrap(~rep, scales="fixed",ncol=6)+		
	theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),
		axis.title.y=element_blank(),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))

###   Section 2 when the delta increase in # invaded per plot leads to delta decrease in natives  ----------
library(dplyr);library(ggpubr)
#detach(package:raster)
load("data analysis2/plots.summary.rda")
load("data analysis2/fia.abundance.ena.rda")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=allplots%>%filter(is.na(degreeOfEstablishment)|degreeOfEstablishment%in%status[1:2])
vars=c("HillEven.native","PD.native","M_Ap.native","FAD.native","M_At.native")
inv=allplots%>%filter(invaded==1)%>%select(pltID,PREV_PLT_CN,PLT_CN,pltID.no.year)%>%distinct
inv.pre=allplots%>%filter(PLT_CN%in%inv$PREV_PLT_CN)%>%select(pltID,PREV_PLT_CN,PLT_CN,invaded)%>%distinct

inv.div=allplots%>%filter(pltID%in%inv$pltID)%>%group_by(pltID)%>%
	summarize(pltID.no.year=unique(pltID.no.year),PLT_CN=unique(PLT_CN),PREV_PLT_CN=unique(PREV_PLT_CN),SR.invasive=length(unique(Accepted_name[which(degreeOfEstablishment%in%status[1:2])])),SR.native=length(unique(Accepted_name[which(is.na(degreeOfEstablishment))])),
	BA.invasive=sum(Abund_weight[which(degreeOfEstablishment%in%status[1:2])]),BA.native=sum(Abund_weight[which(is.na(degreeOfEstablishment))]),
	Biomass.invasive=sum(Biomass[which(degreeOfEstablishment%in%status[1:2])]),Biomass.native=sum(Biomass[which(is.na(degreeOfEstablishment))]))%>%
	left_join(plots.div[,c("pltID",vars)],by="pltID")
inv.pre.div=allplots%>%filter(pltID%in%inv.pre$pltID)%>%group_by(pltID)%>%
	summarize(pltID.no.year=unique(pltID.no.year),PLT_CN=unique(PLT_CN),PREV_PLT_CN=unique(PREV_PLT_CN),SR.invasive=length(unique(Accepted_name[which(degreeOfEstablishment%in%status[1:2])])),SR.native=length(unique(Accepted_name[which(is.na(degreeOfEstablishment))])),
	BA.invasive=sum(Abund_weight[which(degreeOfEstablishment%in%status[1:2])]),BA.native=sum(Abund_weight[which(is.na(degreeOfEstablishment))]),
	Biomass.invasive=sum(Biomass[which(degreeOfEstablishment%in%status[1:2])]),Biomass.native=sum(Biomass[which(is.na(degreeOfEstablishment))]))%>%
	left_join(plots.div[,c("pltID",vars)],by="pltID")
inv.div2=inv.div%>%filter(PREV_PLT_CN%in%inv.pre.div$PLT_CN)
inv.div2[is.na(inv.div2)]=0;inv.pre.div[is.na(inv.pre.div)]=0
delta.pre=do.call(cbind,Map('-', inv.div2[,-c(1:4)], inv.pre.div[,-c(1:4)]))%>%cbind(inv.div2[,'pltID'],.)%>%as.data.frame%>%filter(SR.invasive>0)

library(ggcorrplot)	
#the correlation between each pair of variables is computed, missing values are handled by casewise deletion
corr <- cor(delta.pre[,-1],use="pairwise.complete.obs")[c("SR.invasive","BA.invasive","Biomass.invasive"),c("SR.native","BA.native","Biomass.native",vars)]
p.mat <- cor_pmat(delta.pre[,-1],na.action="na.omit")[c("SR.invasive","BA.invasive","Biomass.invasive"),c("SR.native","BA.native","Biomass.native",vars)]
p1=ggcorrplot(t(corr),method = "square", hc.order = FALSE, type = "full",tl.srt=0, p.mat = t(p.mat),digits=2,
	outline.color = "gray",insig = "pch",lab = TRUE,colors = c("#6D9EC1", "white", "#E46726"),sig.level = 0.05)+
	theme(axis.text = element_text(size=12,color='black'),axis.text.x = element_text(angle=60),legend.position="right")
#map the Plots
library(sf)
ecoregion=st_read("S_USA.EcoMapProvinces/S_USA.EcoMapProvinces.shp")%>% st_transform("+proj=longlat +datum=WGS84")
theme=theme_bw()+
		theme(axis.text = element_text(size=12,color='black'),			
			axis.title = element_blank(),
			panel.background = element_rect(fill = '#619CFF', colour = 'black'))	
plots.sp=allplots%>%filter(pltID%in%delta.pre$pltID)%>% dplyr::select(pltID.no.year,LON,LAT)%>% distinct%>%
	st_as_sf(coords=c("LON","LAT"))%>% st_set_crs("+proj=longlat +datum=WGS84")
p2=ggplot() + geom_sf(data=ecoregion,color="lightgray",fill ="white") +xlim(-100,-67)+ 
		geom_sf(data=plots.sp,color="#F8766D",size=1,show.legend=F)+
		annotate("text",x=-73 , y=25,size=6,
			label = paste0("Plots = ",nrow(plots.sp)),color = "black")+theme		
plt1=ggarrange(p1,p2,widths=c(1,0.5),ncol=2,labels="auto")
	
#scannter plot with quantile regression
# function for labelling; ref:https://stackoverflow.com/questions/65695409/is-there-a-neat-approach-to-label-a-ggplot-plot-with-the-equation-and-other-stat
stat_rq_eqn <- function(formula = y ~ x, tau = 0.5, colour = "red", label.y = 0.9, ...) {
  require(ggpmisc);require(quantreg)
  stat_fit_tidy(method = "rq",
                method.args = list(formula = formula, tau = tau), 
                tidy.args = list(se.type = "nid"),
                mapping = aes(label = sprintf('italic(tau)~"="~%.3f~"; "~italic(P)~"="~%.3g',
                                              after_stat(x_tau),
                                              after_stat(ifelse(x_p.value<0.006,0.01,round(x_p.value,3))))),
                parse = TRUE,
                colour = colour,
                label.y = label.y,
                ...)
}
theme=theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'))		
p2.1=ggplot(data=delta.pre,aes(x=SR.invasive,y=SR.native)) + 
	geom_point(data=delta.pre[delta.pre[,"SR.native"]>0,],shape=21,color="black",fill ="#6D9EC1",size=1.5,alpha=0.5,position=position_jitter(h=0.15,w=0.01))+
	geom_point(data=delta.pre[delta.pre[,"SR.native"]<0,],shape=21,color="black",fill ="#E46726",size=1.5,alpha=0.5,position=position_jitter(h=0.15,w=0.01))+
	geom_point(data=delta.pre[delta.pre[,"SR.native"]==0,],shape=21,color="black",fill ="gray",size=1.5,alpha=0.5,position=position_jitter(h=0.15,w=0.01))+
	geom_hline(yintercept=0,col="darkgray",size=1,linetype="longdash")+
	geom_quantile(aes(linetype = as.factor(..quantile..)),quantiles = c(0.1,0.3,0.5,0.9),linewidth=1,color="red",alpha=0.8,show.legend=F)+
	stat_rq_eqn(tau = 0.1, colour = "black", label.x = 0.8, label.y = 0.15)+
	stat_rq_eqn(tau = 0.3, colour = "black", label.x = 0.8,label.y = 0.1)+
	stat_rq_eqn(tau = 0.5, colour = "black", label.x = 0.8,label.y = 0.05)+
	stat_rq_eqn(tau = 0.9, colour = "black", label.x = 0.8,label.y = 0.0)+
	theme+labs(x="Increased in invasive richness",y="Changes in naive richness")
p2.2=ggplot(data=delta.pre,aes(x=BA.invasive,y=BA.native)) + 
	geom_point(data=delta.pre[delta.pre[,"BA.native"]>0,],shape=21,color="black",fill ="#6D9EC1",size=1.5,alpha=0.5,position=position_jitter(h=0.15,w=0.01))+
	geom_point(data=delta.pre[delta.pre[,"BA.native"]<0,],shape=21,color="black",fill ="#E46726",size=1.5,alpha=0.5,position=position_jitter(h=0.15,w=0.01))+
	geom_point(data=delta.pre[delta.pre[,"BA.native"]==0,],shape=21,color="black",fill ="gray",size=1.5,alpha=0.5,position=position_jitter(h=0.15,w=0.01))+
	geom_hline(yintercept=0,col="darkgray",size=1,linetype="longdash")+
	geom_quantile(aes(linetype = as.factor(..quantile..)),quantiles = c(0.1,0.3,0.5,0.9),linewidth=1,color="red",alpha=0.8,show.legend=F)+
	stat_rq_eqn(tau = 0.1, colour = "black", label.x = 0.8, label.y = 0.15)+
	stat_rq_eqn(tau = 0.3, colour = "black", label.x = 0.8,label.y = 0.1)+
	stat_rq_eqn(tau = 0.5, colour = "black", label.x = 0.8,label.y = 0.05)+
	stat_rq_eqn(tau = 0.9, colour = "black", label.x = 0.8,label.y = 0.0)+
	theme+labs(x="Increased in invasive basal area",y="Changes in naive basal area")
p2.3=ggplot(data=delta.pre,aes(x=Biomass.invasive,y=Biomass.native)) + 
	geom_point(data=delta.pre[delta.pre[,"Biomass.native"]>0,],shape=21,color="black",fill ="#6D9EC1",size=1.5,alpha=0.5,position=position_jitter(h=0.15,w=0.01))+
	geom_point(data=delta.pre[delta.pre[,"Biomass.native"]<0,],shape=21,color="black",fill ="#E46726",size=1.5,alpha=0.5,position=position_jitter(h=0.15,w=0.01))+
	geom_point(data=delta.pre[delta.pre[,"Biomass.native"]==0,],shape=21,color="black",fill ="gray",size=1.5,alpha=0.5,position=position_jitter(h=0.15,w=0.01))+
	geom_hline(yintercept=0,col="darkgray",size=1,linetype="longdash")+
	geom_quantile(aes(linetype = as.factor(..quantile..)),quantiles = c(0.1,0.3,0.5,0.9),linewidth=1,color="red",alpha=0.8)+
	stat_rq_eqn(tau = 0.1, colour = "black", label.x = 0.8, label.y = 0.15)+
	stat_rq_eqn(tau = 0.3, colour = "black", label.x = 0.8,label.y = 0.1)+
	stat_rq_eqn(tau = 0.5, colour = "black", label.x = 0.8,label.y = 0.05)+
	stat_rq_eqn(tau = 0.9, colour = "black", label.x = 0.8,label.y = 0.0)+
	theme+labs(x="Increased in invasive biomass",y="Changes in naive biomass",linetype ="Quantile")+
	theme(legend.position=c(0.7,0.4),legend.background = element_rect(fill = NA),
		legend.title=element_text(face="bold",size=12),
		legend.text=element_text(face="bold",size=10))
plt2=ggarrange(p2.1,p2.2,p2.3,nrow=1,labels=c("c","d","e"))
ggarrange(plt1,plt2,nrow=2)		

#impacts of each invasive species
se=function(x) sd(x, na.rm = TRUE)/length(na.omit(x))
inv.sp=allplots%>%filter(MEASYEAR>=2000&MEASYEAR<=2021)%>%filter(!is.na(degreeOfEstablishment))%>%select(pltID,Accepted_name)%>%distinct
delta.pre2=delta.pre;#delta.pre2[,-c(1:2)]=scale(delta.pre2[,-c(1:2)])
inv.se=delta.pre2%>%left_join(inv.sp,by="pltID")%>%filter(!is.na(Accepted_name))%>%group_by(Accepted_name)%>%
	summarise(across(SR.native:M_At.native, ~ se(.x)))%>%
	tidyr::pivot_longer(cols=SR.native:M_At.native,names_to = "Vars", values_to = "se")
inv.se[is.na(inv.se)]=0	
inv.stat.n=delta.pre2%>%left_join(inv.sp,by="pltID")%>%filter(!is.na(Accepted_name))%>%group_by(Accepted_name)%>%tally
inv.stat=delta.pre2%>%left_join(inv.sp,by="pltID")%>%filter(!is.na(Accepted_name))%>%group_by(Accepted_name)%>%
	summarise(across(SR.native:M_At.native, ~ mean(.x[SR.invasive>0], na.rm = TRUE)))%>%
	tidyr::pivot_longer(cols=SR.native:M_At.native,names_to = "Vars", values_to = "mean")%>%
	left_join(inv.se,by=c("Accepted_name","Vars"))	%>%
	filter(Vars%in%c("SR.native","BA.native","Biomass.native",vars))%>%
	left_join(inv.stat.n,by="Accepted_name")	
inv.stat$Vars=factor(inv.stat$Vars,levels=rev(c("SR.native","BA.native","Biomass.native",vars)))	
inv.stat$type=ifelse((inv.stat$mean+1.96*inv.stat$se)<0,"negative",ifelse((inv.stat$mean-1.96*inv.stat$se)>0,"positive","insig"))
#24 out of 28 species are found here
#4 spp (Populus alba, Tamarindus indica, Salix x fragilis, Leucaena leucocephala) only show up once, 1 sp (Elaeagnus angustifolia) occur in plots without remeasurements
ggplot() +geom_errorbar(data=inv.stat,aes(x=mean, y=Vars,xmin=mean-1.96*se, xmax=mean+1.96*se,color=type),width=1,size=1,show.legend=FALSE,alpha=0.8)+
	geom_point(data=inv.stat,aes(x=mean, y=Vars,fill=type),size=2,shape=21,color="black",show.legend=F) + 		
	scale_fill_manual(values=c("insig"="gray", "positive"="#6D9EC1","negative"="#E46726"))+	
	scale_color_manual(values=c("insig"="gray", "positive"="#6D9EC1","negative"="#E46726"))+	
	labs(x="Difference after invaded") +
	geom_text(data = inv.stat.n, mapping = aes(x=0 , y=-0.5, label = paste0("nPlots = ",n)), hjust   = "center", vjust   = "inward",size=4)+
	facet_wrap(~Accepted_name, scales="fixed",ncol=6)+		
	theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),
		axis.title.y=element_blank(),
		#axis.text.y = element_text(size=6,face="italic"),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))+
	geom_vline(xintercept=0,col="darkgray",size=1,linetype="longdash")

#map plots with inc/dec richness/ba/biomass
theme=theme_bw()+
		theme(axis.text = element_blank(),			
			axis.title = element_blank(),
			panel.background = element_rect(fill = '#619CFF', colour = 'black'),
			legend.position="right",legend.background = element_rect(fill = NA),
		legend.title=element_text(face="bold",size=12),
		legend.text=element_text(face="bold",size=10))	
plots.vaule=delta.pre%>%left_join(allplots[,c("pltID","LON","LAT")],by="pltID")%>%
	st_as_sf(coords=c("LON","LAT"))%>% st_set_crs("+proj=longlat +datum=WGS84")
p1=ggplot() + geom_sf(data=ecoregion,color="darkgray",fill ="lightgray") +xlim(-100,-67)+ 
		geom_sf(aes(color=SR.native),data=plots.vaule,size=1.5,show.legend=T,alpha=0.1,na.rm = T)+
		scale_colour_gradient2(low = "#E46726",  mid = "white",  high = "#6D9EC1",  midpoint = 0)+
		theme+labs(color="Changes in\nNative\nRichness")	
p2=ggplot() + geom_sf(data=ecoregion,color="darkgray",fill ="lightgray") +xlim(-100,-67)+ 
		geom_sf(aes(color=BA.native),data=plots.vaule,size=1.5,show.legend=T,alpha=0.1,na.rm = T)+
		scale_colour_gradient2(low = "#E46726",  mid = "white",  high = "#6D9EC1",  midpoint = 0)+
		theme+labs(color="Changes in\nNative\nBasal Area")	
p3=ggplot() + geom_sf(data=ecoregion,color="darkgray",fill ="lightgray") +xlim(-100,-67)+ 
		geom_sf(aes(color=Biomass.native),data=plots.vaule,size=1.5,show.legend=T,alpha=0.1,na.rm = T)+
		scale_colour_gradient2(low = "#E46726",  mid = "white",  high = "#6D9EC1",  midpoint = 0)+
		theme+labs(color="Changes in\nNative\nBiomass")
p4=ggplot() + geom_sf(data=ecoregion,color="darkgray",fill ="lightgray") +xlim(-100,-67)+ 
		geom_sf(aes(color=HillEven.native),data=plots.vaule,size=1.5,show.legend=T,alpha=0.1,na.rm = T)+
		scale_colour_gradient2(low = "#E46726",  mid = "white",  high = "#6D9EC1",  midpoint = 0)+
		theme+labs(color="Changes in\nNative\nHillEvenness")	
p5=ggplot() + geom_sf(data=ecoregion,color="darkgray",fill ="lightgray") +xlim(-100,-67)+ 
		geom_sf(aes(color=PD.native),data=plots.vaule,size=1.5,show.legend=T,alpha=0.1,na.rm = T)+
		scale_colour_gradient2(low = "#E46726",  mid = "white",  high = "#6D9EC1",  midpoint = 0)+
		theme+labs(color="Changes in\nNative\nPD")	
p6=ggplot() + geom_sf(data=ecoregion,color="darkgray",fill ="lightgray") +xlim(-100,-67)+ 
		geom_sf(aes(color=M_Ap.native),data=plots.vaule,size=1.5,show.legend=T,alpha=0.1,na.rm = T)+
		scale_colour_gradient2(low = "#E46726",  mid = "white",  high = "#6D9EC1",  midpoint = 0)+
		theme+labs(color="Changes in\nNative\nM_Ap")
p7=ggplot() + geom_sf(data=ecoregion,color="darkgray",fill ="lightgray") +xlim(-100,-67)+ 
		geom_sf(aes(color=FAD.native),data=plots.vaule,size=1.5,show.legend=T,alpha=0.1,na.rm = T)+
		scale_colour_gradient2(low = "#E46726",  mid = "white",  high = "#6D9EC1",  midpoint = 0)+
		theme+labs(color="Changes in\nNative\nFAD")	
p8=ggplot() + geom_sf(data=ecoregion,color="darkgray",fill ="lightgray") +xlim(-100,-67)+ 
		geom_sf(aes(color=M_At.native),data=plots.vaule,size=1.5,show.legend=T,alpha=0.1,na.rm = T)+
		scale_colour_gradient2(low = "#E46726",  mid = "white",  high = "#6D9EC1",  midpoint = 0)+
		theme+labs(color="Changes in\nNative\nM_At")
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,labels="auto")		


##compare diversity and biomass between invaded and the uninvaded plots within each ecoregion
## the invaded plots are plots that are invaded in future measurements while not being invaded yet
library(dplyr);library(ggpubr)
load("data analysis2/plots.summary.rda")
load("data analysis2/fia.abundance.rda")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=fia%>%filter((!is.na(Abund_weight))&LON>=-100&is.na(degreeOfEstablishment)|degreeOfEstablishment%in%status[1:2])%>%group_by(pltID)%>%
	summarize(SR.native=length(unique(Accepted_name[which(is.na(degreeOfEstablishment))])),BA.native=sum(Abund_weight[which(is.na(degreeOfEstablishment))]),
			Biomass.native=sum(Biomass[which(is.na(degreeOfEstablishment))]),
			SR.invasive=length(unique(Accepted_name[which(!is.na(degreeOfEstablishment))])),BA.invasive=sum(Abund_weight[which(!is.na(degreeOfEstablishment))]),
			Biomass.invasive=sum(Biomass[which(!is.na(degreeOfEstablishment))]),LON=unique(LON),LAT=unique(LAT))
vars=c("HillEven.native","PD.native","M_Ap.native","FAD.native","M_At.native")
plots.div=plots.div%>%select(year,ecoregion,pltID,pltID.no.year,invaded,matches(vars))%>%left_join(allplots,by="pltID")%>%filter(!is.na(SR.native))

inv.plots=plots.div%>%filter(invaded==1)%>%select(pltID.no.year)%>%distinct
inv.plots$inasion.future=1
plots.div=plots.div%>%left_join(inv.plots,by="pltID.no.year")
plots.div[is.na(plots.div$inasion.future),"inasion.future"]=0
plots.div$inv.type=ifelse(plots.div$inasion.future==0,"Native",ifelse(plots.div$invaded==0,"To be invaded","Invaded"))
plots.div%>%group_by(invaded,inasion.future,inv.type)%>%tally#summarise(n=length(unique(pltID.no.year)))
#plots.div[,c(vars,"SR.native","BA.native","Biomass.native")]=scale(plots.div[,c(vars,"SR.native","BA.native","Biomass.native")])

se=function(x) sd(x, na.rm = TRUE)/length(na.omit(x))
stat.se=plots.div%>%group_by(ecoregion,inv.type)%>%
	summarise(across(HillEven.native:Biomass.invasive, ~ se(.x)),.groups="keep")%>%
	tidyr::pivot_longer(cols=HillEven.native:Biomass.invasive,names_to = "Vars", values_to = "se")
stat.all=plots.div%>%group_by(ecoregion,inv.type)%>%
	summarise(nPlots=length(unique(pltID.no.year)),across(HillEven.native:Biomass.invasive, ~ mean(.x,na.rm=T)),.groups="keep")%>%
	tidyr::pivot_longer(cols=HillEven.native:Biomass.invasive,names_to = "Vars", values_to = "mean")%>%
	left_join(stat.se,by=c("ecoregion","inv.type","Vars"))
stat.all[is.na(stat.all)]=0
low.eco=stat.all%>%filter((inv.type=="To be invaded"&nPlots<5)|(inv.type=="Native"&nPlots<5))%>%select(ecoregion)%>%distinct
stat.all=stat.all%>%filter(!ecoregion%in%low.eco$ecoregion)
inv.stat.n=stat.all%>%select(ecoregion,inv.type,nPlots)%>%distinct%>%tidyr::pivot_wider(names_from = inv.type, values_from = nPlots)%>% tibble::column_to_rownames("ecoregion")
write.csv(inv.stat.n,"data analysis2/inv.stat.csv")

stat=stat.all%>%filter(Vars%in%c("SR.native","Biomass.native","M_Ap.native"))
stat$Vars=factor(stat$Vars,levels=c("SR.native","Biomass.native","M_Ap.native"))
ggplot() +geom_errorbar(data=stat,aes(x=mean, y=ecoregion,xmin=mean-1.96*se, xmax=mean+1.96*se,color=inv.type),width=0.8,size=1,show.legend=FALSE,alpha=0.5)+
	geom_point(data=stat,aes(x=mean, y=ecoregion,fill=inv.type),size=1.4,shape=21,color="black",show.legend=T,alpha=0.8) + 		
	labs(x="Mean values",y="Ecoprovinces") +
	guides(color=guide_legend(ncol=1,byrow=T,override.aes = list(size = 4)))+	
	# geom_text(data = inv.stat.n%>%filter(inv.type=="Native"), mapping = aes(x=0 , y=-0.5,color=inv.type, label = paste0("nPlots = ",nPlots)), hjust   = "center", vjust   = "inward",size=4,show.legend=F)+
	# geom_text(data = inv.stat.n%>%filter(inv.type=="To be invaded"), mapping = aes(x=0 , y=-1,color=inv.type, label = paste0("nPlots = ",nPlots)), hjust   = "center", vjust   = "inward",size=4,show.legend=F)+
	scale_fill_manual(values=c("Native"="#619CFF", "To be invaded"="#F8766D","Invaded"="darkred"))+
	scale_color_manual(values=c("Native"="#619CFF", "To be invaded"="#F8766D","Invaded"="darkred"))+
	facet_wrap(~Vars, scales="free_x")+		
	theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),
		legend.position = "top",
		#axis.text.y = element_text(size=6,face="italic"),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))

library(sf);
#map the Plots
ecoregion=st_read("S_USA.EcoMapProvinces/S_USA.EcoMapProvinces.shp")%>% st_transform("+proj=longlat +datum=WGS84")
the=theme_bw()+theme(panel.background = element_rect(fill = NA, colour = 'black'),
			axis.text = element_text(size=12,color='black'),
			axis.title = element_blank(),
			panel.grid =element_blank(),
			#legend.position = c(0.85,0.2),
			legend.background=element_rect(fill='transparent'),
			legend.text=element_text(face="bold",size=10),
			legend.title=element_text(face="bold",size=10),
			legend.key = element_blank())
ggplot() +  
		geom_sf(data=fia.points,aes(color=ecoregion),size=1)+
		geom_sf(data=ecoregion,color="lightgray",fill ="transparent") +xlim(-100,-67)+
		scale_color_manual(values=c("#00BA38","green","brown","darkgreen","#93AA00","#00C19F","#F8766D","#DB72FB","#FF61C3","darkgray",
		"black","purple","orange","lightgreen","pink"))+
		guides(color=guide_legend(ncol=1,byrow=T,override.aes = list(size = 3)))+	
		the	
		
#Following Jeremy's suggestion on Thu 5/9/2024		
# defines the data subset is the appearance of exactly one individual belonging to a species that was not present at t1: Identify all FIA plots where one individual of a "new" species appeared between t1 and t2. 
# after controlling for other covariates, does it matter if the new species is native or non-native? 
# Specifically, quantify how the response variables (changes in diversity metrics and biomass from t2 to t3) depend on explanatory variables: 
# stand age at t1, biomass at t1, change in biomass from t1 to t2, biomass lost due to harvest from t2 to t3, and whether the new species is native or non-native. 	
# To control for climate and geography, separate multiple regression models for each ecoprovince	
library(dplyr)
load("data analysis2/fia.abundance.rda")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=fia%>%filter((!is.na(Abund_weight))&LON>=-100&(is.na(degreeOfEstablishment)|degreeOfEstablishment%in%status[1:2]))
allplots$PREV_TRE_CN2=ifelse(is.na(allplots$PREV_TRE_CN),NA,"Not NA")
plots=allplots%>%filter(!is.na(PREV_PLT_CN))#%>%filter(MEASYEAR>=2000&MEASYEAR<=2021)

splist=sort(unique(plots$Accepted_name))
sp.status=plots%>%select(Accepted_name,degreeOfEstablishment)%>%distinct
sp.status$invasive=ifelse(is.na(sp.status$degreeOfEstablishment),0,1)
inv.plots=allplots%>%filter(degreeOfEstablishment%in%status[1:2])%>%select(PLT_CN)%>%distinct()
delta.change=c()	
for (sp in splist){
	plotID.t2=allplots%>%filter(Accepted_name==sp&(!is.na(PREV_PLT_CN)))%>%select(PLT_CN,PREV_PLT_CN)%>%distinct()
	plotID.t1=allplots%>%filter((Accepted_name!=sp)&(PREV_PLT_CN%in%plotID.t2$PREV_PLT_CN))%>%select(PLT_CN)%>%distinct()
	#plotID.t3=allplots%>%filter((Accepted_name==sp)&(PREV_PLT_CN%in%plotID.t2$PLT_CN))%>%select(PREV_PLT_CN)%>%distinct()
	plotID.t2=plotID.t2%>%filter(PREV_PLT_CN%in%plotID.t1$PLT_CN&PLT_CN%in%plotID.t3$PREV_PLT_CN)
	invasive=sp.status[sp.status$Accepted_name==sp,"invasive"]
	#if (invasive==0) plotID.t2=plotID.t2%>%filter((!PLT_CN%in%inv.plots$PLT_CN)&(!PLT_CN%in%inv.plots$PLT_CN))
	if(nrow(plotID.t2)==0) next;
	
	plots.t1=allplots%>%filter(PLT_CN%in%plotID.t2$PREV_PLT_CN)%>%group_by(PLT_CN)%>%
		summarize(STDAGE=mean(STDAGE,na.rm=T),ECOSUBCD=unique(ECOSUBCD),Biomass.t1=sum(Biomass),
			SR.t1=length(unique(Accepted_name)),BA.t1=sum(Abund_weight))		
	plots.t2=allplots%>%filter(PLT_CN%in%plotID.t2$PLT_CN)%>%group_by(PLT_CN,PREV_PLT_CN)%>%
		summarize(Biomass.t2=sum(Biomass[which(Accepted_name!=sp)]),SR.t2=length(unique(Accepted_name))-1,BA.t2=sum(Abund_weight[which(Accepted_name!=sp)]),.groups="keep")	
	plots.t3=allplots%>%filter(PREV_PLT_CN%in%plots.t2$PLT_CN)%>%group_by(PREV_PLT_CN)%>%
		summarize(Biomass.t3=sum(Biomass[which(Accepted_name!=sp)]),SR.t3=length(unique(Accepted_name))-1,BA.t3=sum(Abund_weight[which(Accepted_name!=sp)]))
	plots.t=plots.t1%>%left_join(plots.t2,by=c("PLT_CN"="PREV_PLT_CN"))	%>%left_join(plots.t3,by=c("PLT_CN.y"="PREV_PLT_CN"))%>%
		mutate(SR.t12.por=(SR.t1-SR.t2)/SR.t1, BA.t12.por=(BA.t1-BA.t2)/BA.t1, Biomass.t12.por=(Biomass.t1-Biomass.t2)/Biomass.t1,
			SR.t23.por=(SR.t2-SR.t3)/SR.t2, BA.t23.por=(BA.t2-BA.t3)/BA.t2, Biomass.t23.por=(Biomass.t2-Biomass.t3)/Biomass.t2,
			SR.t13.por=(SR.t1-SR.t3)/SR.t1, BA.t13.por=(BA.t1-BA.t3)/BA.t1, Biomass.t13.por=(Biomass.t1-Biomass.t3)/Biomass.t1,
			
			SR.t12=SR.t1-SR.t2, BA.t12=BA.t1-BA.t2, Biomass.t12=Biomass.t1-Biomass.t2,
			SR.t23=SR.t2-SR.t3, BA.t23=BA.t2-BA.t3, Biomass.t23=Biomass.t2-Biomass.t3,
			SR.t13=SR.t1-SR.t3, BA.t13=BA.t1-BA.t3, Biomass.t13=Biomass.t1-Biomass.t3)	
	plots.t$Accepted_name=sp;	plots.t$invasive=invasive;	plots.t$nplots=nrow(plotID.t2)
	delta.change=rbind(delta.change,plots.t)	
}
delta.change$invasive=factor(delta.change$invasive,levels=c(0,1))
plot(delta.change[,c("invasive","SR.t23")])

library(lme4);
#library(afex)#This will automatically add a p-value column to the output of the lmer for the fixed effects.
library(ggplot2)

eco=unique(delta.change$ECOSUBCD)
stat.all=stat.n.all=c()
for (i in eco){
	dat=delta.change%>%filter(ECOSUBCD%in%i)%>%select(Accepted_name,invasive,nplots,STDAGE,SR.t1,BA.t1,Biomass.t1,Biomass.t12,SR.t23)%>%
			na.omit%>%filter_all(all_vars(!is.infinite(.)))%>%filter(nplots>10)		
	stat.n=dat%>%group_by(invasive)%>%tally
	if (nrow(stat.n)==1) next;
	if (stat.n[stat.n$invasive==1,"n"]<10) next; 
	dat[,c("STDAGE","SR.t1","BA.t1","Biomass.t1","Biomass.t12")]=scale(dat[,c("STDAGE","SR.t1","BA.t1","Biomass.t1","Biomass.t12")])
	M <- lmer(SR.t23 ~ invasive + STDAGE + SR.t1 + BA.t1 + Biomass.t1 + Biomass.t12 +(1|Accepted_name) ,data=dat)
	stat=as.data.frame(summary(M)$coef[-1,c(1,2,5)])
	colnames(stat)[2:3]=c("se","p")
	stat$var=rownames(stat)
	stat$ecoregion=stat.n$ecoregion=i
	stat.all=rbind(stat.all,stat)
	stat.n.all=rbind(stat.n.all,stat.n)	
}
stat.all$type=ifelse((stat.all$Estimate+1.96*stat.all$se)<0,"negative",ifelse((stat.all$Estimate-1.96*stat.all$se)>0,"positive","insig"))
stat.all$var=factor(stat.all$var,levels=rev(c("invasive1","SR.t1","BA.t1","Biomass.t1","Biomass.t12","STDAGE")))
ggplot() +geom_errorbar(data=stat.all,aes(x=Estimate, y=var,xmin=Estimate-1.96*se, xmax=Estimate+1.96*se,color=type),width=1,size=1,show.legend=FALSE,alpha=0.8)+
	geom_point(data=stat.all,aes(x=Estimate, y=var,fill=type),size=2,shape=21,color="black",show.legend=F) + 		
	scale_fill_manual(values=c("insig"="darkgray", "positive"="#6D9EC1","negative"="#E46726"))+	
	scale_color_manual(values=c("insig"="darkgray", "positive"="#6D9EC1","negative"="#E46726"))+	
	labs(x="Estimates") +
	geom_text(data = stat.n.all%>%filter(invasive==0)%>%select(ecoregion,n), mapping = aes(x=0 , y=-0.5, label = paste0("Plots (Native) = ",n)), hjust   = "center", vjust   = "inward",size=4)+
	geom_text(data = stat.n.all%>%filter(invasive==1)%>%select(ecoregion,n), mapping = aes(x=0 , y=-1.2, label = paste0("Plots (Invasive) = ",n)), hjust   = "center", vjust   = "inward",size=4)+	
	facet_wrap(~ecoregion, scales="fixed")+		
	theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),
		axis.title.y=element_blank(),
		#axis.text.y = element_text(size=6,face="italic"),
		strip.text = element_text(colour = 'black', face="italic",size = rel(1.2)), 
		strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))+
	geom_vline(xintercept=0,col="lightgray",size=0.8,linetype="longdash")
#eco2=unique(stat.all$ecoregion)
#for all plots
dat=delta.change%>%select(Accepted_name,ECOSUBCD,invasive,nplots,STDAGE,SR.t1,BA.t1,Biomass.t1,Biomass.t12,SR.t23)%>%
			na.omit%>%filter_all(all_vars(!is.infinite(.)))%>%filter(nplots>10)#%>%filter(nplots>10&ECOSUBCD%in%eco2)		
stat.n=dat%>%group_by(invasive)%>%tally
dat[,c("STDAGE","SR.t1","BA.t1","Biomass.t1","Biomass.t12")]=scale(dat[,c("STDAGE","SR.t1","BA.t1","Biomass.t1","Biomass.t12")])
#dat$Biomass.t23=-dat$Biomass.t23
M <- lmer(SR.t23 ~ invasive + STDAGE + SR.t1 + BA.t1 + Biomass.t1 + Biomass.t12 +(1|Accepted_name)+(1|ECOSUBCD) ,data=dat)
stat=as.data.frame(summary(M)$coef[-1,c(1,2,5)])
colnames(stat)[2:3]=c("se","p")
stat$var=rownames(stat)
stat$type=ifelse((stat$Estimate+1.96*stat$se)<0,"negative",ifelse((stat$Estimate-1.96*stat$se)>0,"positive","insig"))
stat$var=factor(stat$var,levels=rev(c("invasive1","SR.t1","BA.t1","Biomass.t1","Biomass.t12","STDAGE")))

ggplot() +geom_errorbar(data=stat,aes(x=Estimate, y=var,xmin=Estimate-1.96*se, xmax=Estimate+1.96*se,color=type),width=1,size=1,show.legend=FALSE,alpha=0.8)+
	geom_point(data=stat,aes(x=Estimate, y=var,fill=type),size=2,shape=21,color="black",show.legend=F) + 		
	scale_fill_manual(values=c("insig"="darkgray", "positive"="#6D9EC1","negative"="#E46726"))+	
	scale_color_manual(values=c("insig"="darkgray", "positive"="#6D9EC1","negative"="#E46726"))+	
	labs(x="Estimates") +
	geom_text(data = stat.n%>%filter(invasive==0), mapping = aes(x=0 , y=-0.5, label = paste0("Plots (Native) = ",n)), hjust   = "center", vjust   = "inward",size=4)+
	geom_text(data = stat.n%>%filter(invasive==1), mapping = aes(x=0 , y=-1.2, label = paste0("Plots (Invasive) = ",n)), hjust   = "center", vjust   = "inward",size=4)+	
	theme_bw()+theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),
		axis.title.y=element_blank())+
	geom_vline(xintercept=0,col="lightgray",size=0.8,linetype="longdash")

##count how many times different non-native situations occur
##Specifically, how many plot locations do the following situations occur
# One (and only 1) non-native species appears between t1 and t2, with only 1 individual present at t2.
# One (and only 1) non-native species appears between t1 and t2, with > 1 individual present at t2.
# More than one non-native species appear between t1 and t2. 

# You could further divide the cases based on the following:
# Zero non-native species present at t1.
# At least one non-native species present at t1. 

# There are 3 x 2 = 6 combinations of the above. 
# within case 2b, it would be interesting to make a table of whether the (only) new individual of the invasive species appears as a small tree in the microplot or a large tree in the subplot
# All new stems will be either RECONCILECD = 1  for ingrowth or RECONCILECD = 2  for through growth.  (See definitions from FIA user guide below)
# If the tree is RECONCILECD = 1 and the DBH is less than 5 inches / 12.5 cm, we know it is a new sapling in the microplot.
# If the tree is RECONCILECD = 1 and the DBH is greater than 5 inches / 12.5 cm we know it is a new tree, coming into the census in the subplot  
# For trees that go from <5 inches / 12.5 to  greater than that size in the microplot, those will have RECONCILECD = 2  “through growth”
library(dplyr)
load("data analysis2/fia.abundance.rda")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=fia%>%filter(STATUSCD==1&LON>=-100&CONDPROP_UNADJ==1)
#a=allplots%>%filter(!is.na(PREV_PLT_CN))%>%plots%>%group_by(RECONCILECD,DIA2,MICR_COMPONENT_AL_FOREST,SUBP_COMPONENT_AL_FOREST,PREV_TRE_CN2)%>%tally
#write(a,"ingrowth.csv")
inv.sp=allplots%>%filter(degreeOfEstablishment%in%status[1:2])%>%select(Accepted_name)%>%distinct
inv.plot=allplots%>%filter(degreeOfEstablishment%in%status[1:2])%>%select(PLT_CN)%>%distinct
re.all=c()
for(plt in inv.plot$PLT_CN){
# plt=inv.plot$PLT_CN[1324]
	t2=allplots%>%filter(PLT_CN%in%plt)
	if (is.na(unique(t2$PREV_PLT_CN))) next;
	t1=allplots%>%filter(PLT_CN%in%unique(t2$PREV_PLT_CN))
	#new.inv.t2=t2%>%filter(RECONCILECD%in%c(1,2)&Accepted_name%in%inv.sp$Accepted_name)
	# new.inv.t2=t2%>%filter((Accepted_name%in%inv.sp$Accepted_name)&
		# (SUBP_COMPONENT_AL_FOREST=='INGROWTH'|MICR_COMPONENT_AL_FOREST=='INGROWTH')&
		# (!SUBP_COMPONENT_AL_FOREST%in%c("REVERSION1","SURVIVOR","REVERSION2"))&(!MICR_COMPONENT_AL_FOREST%in%c("REVERSION1","SURVIVOR","REVERSION2")))
	new.inv.t2=t2%>%filter(is.na(PREV_TRE_CN)&Accepted_name%in%inv.sp$Accepted_name)	
	case.t1=ifelse(sum(t1$Accepted_name%in%inv.sp$Accepted_name)>0,">=1 inv.sp in t1","0 inv.sp in t1")
	case.t2=ifelse(nrow(new.inv.t2)==0,"No new inv.sp in t2",
				ifelse(nrow(new.inv.t2)==1,"1 new inv.sp with 1 new individual in t2",
					ifelse(length(unique(new.inv.t2$Accepted_name))==1,"1 new inv.sp with >1 new individual in t2",">1 new inv.sp in t2")))
	if (nrow(new.inv.t2)>0) {
		mf2=function(x) paste0(x,collapse =" individual ")
		#new.inv.t2$new.tree=ifelse(new.inv.t2$DIA<12.5,"new sapling",ifelse((new.inv.t2$RECONCILECD==2)|(is.na(new.inv.t2$RECONCILECD)),"new tree (grew from saplings in t1)","new tree"))
		new.inv.t2$new.tree=ifelse(new.inv.t2$DIA<12.5,"new sapling","new tree")		
		stat=new.inv.t2%>%group_by(new.tree)%>%tally
		stat$n=ifelse(stat$n==1,"=1",">1")
	}
	case.t2.note=ifelse(nrow(new.inv.t2)==0,NA,paste(apply(stat,1,mf2),collapse=";"))
	re=data.frame(PLT_CN=plt,case.t1=case.t1,case.t2=case.t2,type.of.newtree.t2=case.t2.note)
	re.all=rbind(re.all,re)
	print(which(inv.plot$PLT_CN==plt))
}
write.csv(re.all,"data analysis2/stat.newtree.csv")

stat.all=re.all%>%group_by(case.t1,case.t2,type.of.newtree.t2)%>%tally
stat.all$por=stat.all$n/sum(stat.all$n)*100
write.csv(stat.all,"data analysis2/stats.csv")

##regressions on plots that only one speceis is new in t2 and lasted until t3
#test effect of filters on sample size
test.filter=function(plt,allplots,sp.status){
	require(dplyr)
	t2=allplots%>%filter(PLT_CN%in%plt)
	t3=allplots%>%filter(PREV_PLT_CN%in%unique(t2$PLT_CN))
	case1=case2=case3=case4=NA
	if (nrow(t3)==0) case1="No plot mears or not used in t3"
	new.inv.t2=t2%>%filter(is.na(PREV_TRE_CN))%>%distinct
	if (length(unique(new.inv.t2$Accepted_name))>1) case2=">1 new sp in t2"	
	if (length(unique(new.inv.t2$Accepted_name))==0) case2="No new sp in t2"	
	if (length(unique(new.inv.t2$Accepted_name))==1){
		status=sp.status[sp.status$Accepted_name%in%new.inv.t2$Accepted_name,"degreeOfEstablishment"]
		case3=ifelse(is.na(status),"only 1 new sp and it is native",ifelse(sum(t3$PREV_TRE_CN%in%new.inv.t2$CN)==0,"only 1 new nonative sp in t2 yet no survivor in t3",NA))	
	}
	t1=allplots%>%filter(PLT_CN%in%unique(t2$PREV_PLT_CN))	
	if (nrow(t1)==0) case4="Plot mears in t1 not used"
	re=data.frame(PLT_CN=plt,case1,case2,case3,case4)
	return(re)
}
plots.cn=allplots%>%filter(pltID.no.year%in%plots$pltID.no.year&(degreeOfEstablishment%in%status[1:2]))%>%select(PLT_CN,PREV_PLT_CN)%>%distinct%>%na.omit
re.all=do.call(rbind,parLapply(cl=mycl,plots.cn$PLT_CN,test.filter,allplots,sp.status))
a=re.all%>%group_by(case1,case4,case2,case3)%>%tally
write.csv(a,"a.csv")

#formal analysis
get.plt.reg=function(plt,allplots,sp.status,tree,trait_mat){
#re.all=c()
#for(plt in plots.cn$PLT_CN[which(plots.cn$PLT_CN==plt):nrow(plots.cn)]){
# plt=plots.cn$PLT_CN[1324]
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
	if (nrow(t3)==0) return(NULL) # next;
	new.inv.t2=t2%>%filter(is.na(PREV_TRE_CN))%>%distinct
	if (length(unique(new.inv.t2$Accepted_name))!=1) return(NULL) #next;	
	if (sum(t3$PREV_TRE_CN%in%new.inv.t2$CN)==0) return(NULL) #next;	
	t1=allplots%>%filter(PLT_CN%in%unique(t2$PREV_PLT_CN))	
	if (nrow(t1)==0) return(NULL) #next;
	
	measure.time12=unique(t2$measure.time)-unique(t1$measure.time)
	measure.time23=unique(t3$measure.time)-unique(t2$measure.time)
	new.sp=new.inv.t2$Accepted_name
	case.t1=ifelse(sum(t1$Accepted_name%in%new.sp)==0,"INVADER","SURVIVOR")		
	# if (case.t1=="SURVIVOR"){
		# vars.t1=t1%>%group_by(PLT_CN)%>%
			# summarize(Biomass.t1=sum(Biomass[!Accepted_name%in%new.sp]),SR.t1=length(unique(Accepted_name))-1,BA.t1=sum(Abund_weight[!Accepted_name%in%new.sp]),
				# Biomass.inv.t1=sum(Biomass[Accepted_name%in%new.sp]),BA.inv.t1=sum(Abund_weight[Accepted_name%in%new.sp]))
	# }else{
		# vars.t1=t1%>%group_by(PLT_CN)%>%
			# summarize(Biomass.t1=sum(Biomass[!Accepted_name%in%new.sp]),SR.t1=length(unique(Accepted_name))-1,BA.t1=sum(Abund_weight[!Accepted_name%in%new.sp]),
				# Biomass.inv.t1=0,BA.inv.t1=0)
	# }
	vars.t2=t2%>%group_by(PLT_CN,PREV_PLT_CN)%>%
		summarize(STDAGE=mean(STDAGE,na.rm=T),ECOSUBCD=unique(ECOSUBCD),
			Biomass.t2=sum(Biomass[!Accepted_name%in%new.sp]),SR.t2=length(unique(Accepted_name))-1,BA.t2=sum(Abund_weight[!Accepted_name%in%new.sp]),
			Biomass.inv.t2=sum(Biomass[Accepted_name%in%new.sp]),BA.inv.t2=sum(Abund_weight[Accepted_name%in%new.sp]),.groups="keep")			
	vars.t3=t3%>%group_by(PREV_PLT_CN)%>%
		summarize(Biomass.t3=sum(Biomass[!Accepted_name%in%new.sp]),SR.t3=length(unique(Accepted_name))-1,BA.t3=sum(Abund_weight[!Accepted_name%in%new.sp]),
			Biomass.inv.t3=sum(Biomass[Accepted_name%in%new.sp]),BA.inv.t3=sum(Abund_weight[Accepted_name%in%new.sp]))		
	delta.t23=vars.t2%>%left_join(vars.t3,by=c("PLT_CN"="PREV_PLT_CN"))	%>%
		mutate(SR.t23=(SR.t3-SR.t2)/measure.time23, BA.t23=(BA.t3-BA.t2)/measure.time23, Biomass.t23=(Biomass.t3-Biomass.t2)/measure.time23,
			BA.inv.t23=(BA.inv.t3-BA.inv.t2)/measure.time23, Biomass.inv.t23=(Biomass.inv.t3-Biomass.inv.t2)/measure.time23)	
	
	div.t2 <- t2%>%filter(!Accepted_name%in%new.sp) %>% group_by(PLT_CN,Accepted_name)  %>% summarize(abund_w = length(Abund_weight),.groups ="keep")%>% 
		tidyr::pivot_wider(names_from = Accepted_name, values_from = abund_w) %>% tibble::column_to_rownames("PLT_CN")	%>%div.cal(tree,trait_mat)
	div.t3 <- t3%>%filter(!Accepted_name%in%new.sp) %>% group_by(PREV_PLT_CN,Accepted_name)  %>% summarize(abund_w = length(Abund_weight),.groups ="keep")%>% 
		tidyr::pivot_wider(names_from = Accepted_name, values_from = abund_w) %>% tibble::column_to_rownames("PREV_PLT_CN")	%>%div.cal(tree,trait_mat)
	div.t23=(div.t3-div.t2)/measure.time23
	colnames(div.t23)=paste0(colnames(div.t23),".t23")
	re=data.frame(case.t1=case.t1,Accepted_name=new.sp,degreeOfEstablishment=sp.status[sp.status$Accepted_name%in%new.sp,"degreeOfEstablishment"],
		measure.time12,measure.time23,delta.t23,div.t2,div.t23)
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


## impact of native/non native species on their plots ----
library(dplyr)
load("data analysis2/fia.abundance.rda")
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
allplots=fia%>%filter((!is.na(Abund_weight))&LON>=-100&(is.na(degreeOfEstablishment)|degreeOfEstablishment%in%status[1:2]))
load("data analysis2/plots.summary.rda")
plots=allplots%>%filter(!is.na(PREV_PLT_CN))#%>%filter(MEASYEAR>=2000&MEASYEAR<=2021)
# #check plots that the PREV_PLT_CN have muti PLT_CN
# stat=plots%>%select(PREV_PLT_CN,PLT_CN)%>%distinct%>%group_by(PREV_PLT_CN)%>%tally%>%filter(n>1)
# plots%>%filter(PREV_PLT_CN%in%stat$PREV_PLT_CN)%>%select(pltID,MEASYEAR,pltID.no.year,PLT_CN,PREV_PLT_CN)%>%distinct
splist=unique(plots$Accepted_name)
sp.status=plots%>%select(Accepted_name,degreeOfEstablishment)%>%distinct
sp.status$invasive=ifelse(is.na(sp.status$degreeOfEstablishment),0,1)
vars=c("PD.native","M_Ap.native","FAD.native","M_At.native")
delta.change=c()
for (sp in splist){
	plots.current=allplots%>%filter(Accepted_name==sp&(!is.na(PREV_PLT_CN)))
	plots.pre=allplots%>%filter(PLT_CN%in%plots.current$PREV_PLT_CN)
	if(nrow(plots.pre)==0) next;
	plots.current=allplots%>%filter(PREV_PLT_CN%in%unique(plots.pre$PLT_CN))
	# b=unique(plots.current$PREV_PLT_CN)
	# a=unique(plots.pre$PLT_CN)
	# miss=b[!b%in%a]
	current=plots.current%>%group_by(PREV_PLT_CN)%>%
		summarize(pltID=unique(pltID),SR.native=length(unique(Accepted_name[which(is.na(degreeOfEstablishment))])),BA.native=sum(Abund_weight[which(is.na(degreeOfEstablishment))]),
			Biomass.native=sum(Biomass[which(is.na(degreeOfEstablishment))]))%>%
		left_join(plots.div[,c("pltID",vars)],by="pltID")%>%.[order(.$PREV_PLT_CN),]
	pre=plots.pre%>%group_by(PLT_CN)%>%
		summarize(pltID=unique(pltID),SR.native=length(unique(Accepted_name[which(is.na(degreeOfEstablishment))])),BA.native=sum(Abund_weight[which(is.na(degreeOfEstablishment))]),
			Biomass.native=sum(Biomass[which(is.na(degreeOfEstablishment))]))%>%
		left_join(plots.div[,c("pltID",vars)],by="pltID")%>%.[order(.$PLT_CN),]
	#whether or not the species also occured in previous plots
	plots.pre$inv=ifelse(plots.pre$Accepted_name==sp,1,0)
	invORnot=plots.pre%>%select(PLT_CN,inv)%>%distinct%>%group_by(PLT_CN)%>%summarise(PresentedInPreMeasure=sum(inv))
	#% of changes of biomass and diversity after the plots are invaded by a species
	#mf=function(a,b) (a-b)/b
	delta=do.call(cbind,Map('-', current[,-c(1:2)], pre[,-c(1:2)]))%>%cbind(pre[,'PLT_CN'],.)%>%left_join(invORnot,by="PLT_CN")
	delta.change=rbind(delta.change,data.frame(Accepted_name=sp,delta))
	print(which(splist==sp))
}
delta.change=delta.change%>%left_join(sp.status,by="Accepted_name")
save(delta.change,file="data analysis2/delta.change.rda")

#changes in native diversity and biomass
load("data analysis2/delta.change.rda")
#delta.change[,c(3:9)]=scale(delta.change[,c(3:9)])
#compare based on each species	
get.stat=function(delta.change){
	se=function(x) sd(x, na.rm = TRUE)/length(na.omit(x))
	inv.se=delta.change%>%group_by(Accepted_name)%>% summarise(invasive=unique(invasive),across(SR.native:M_At.native, ~ se(.x)))%>%
		tidyr::pivot_longer(cols=SR.native:M_At.native,names_to = "Vars", values_to = "se")
	inv.se[is.na(inv.se)]=0	
	inv.stat.n=delta.change%>%group_by(Accepted_name)%>%tally
	inv.stat=delta.change%>%group_by(Accepted_name)%>%
		summarise(invasive=unique(invasive),across(SR.native:M_At.native, ~ mean(.x, na.rm = TRUE)))%>%
		tidyr::pivot_longer(cols=SR.native:M_At.native,names_to = "Vars", values_to = "mean")%>%
		left_join(inv.se,by=c("Accepted_name","invasive","Vars"))	%>%
		left_join(inv.stat.n,by="Accepted_name")	# %>% filter(se!=0)
	inv.stat$sig=ifelse((inv.stat$mean+1.96*inv.stat$se)<0,"negative",ifelse((inv.stat$mean-1.96*inv.stat$se)>0,"positive","insig"))
	return(na.omit(inv.stat))
}
stat.n=delta.change%>%filter(PresentedInPreMeasure==0)%>%get.stat
#stat.n%>%group_by(invasive,Vars,sig)%>%tally%>%as.data.frame
stat.y=delta.change%>%filter(PresentedInPreMeasure==1)%>%get.stat

inv.stat%>%na.omit
#compare between invasive and non-invasive
inv.se=delta.change%>%select(-Accepted_name)%>%distinct%>%group_by(invasive)%>% 
	summarise(across(SR.native:M_At.native, ~ se(.x)))%>%
	tidyr::pivot_longer(cols=SR.native:M_At.native,names_to = "Vars", values_to = "se")
inv.se[is.na(inv.se)]=0	
inv.stat=delta.change%>%select(-Accepted_name)%>%distinct%>%group_by(invasive)%>%
	summarise(across(SR.native:M_At.native, ~ mean(.x, na.rm = TRUE)))%>%
	tidyr::pivot_longer(cols=SR.native:M_At.native,names_to = "Vars", values_to = "mean")%>%
	left_join(inv.se,by=c("invasive","Vars"))
## changes after invasion are removed
inv.after=allplots%>%filter(PREV_PLT_CN%in%inv$PLT_CN&invaded==0)%>%select(pltID,PREV_PLT_CN,PLT_CN,invaded)%>%distinct
inv.after.div=allplots%>%filter(pltID%in%inv.after$pltID)%>%group_by(pltID)%>%
	summarize(PLT_CN=unique(PLT_CN),PREV_PLT_CN=unique(PREV_PLT_CN),SR.invasive=length(unique(Accepted_name[which(degreeOfEstablishment%in%status[1:2])])),SR.native=length(unique(Accepted_name[which(is.na(degreeOfEstablishment))])),
	BA.invasive=sum(Abund_weight[which(degreeOfEstablishment%in%status[1:2])]),BA.native=sum(Abund_weight[which(is.na(degreeOfEstablishment))]),
	Biomass.invasive=sum(Biomass[which(degreeOfEstablishment%in%status[1:2])]),Biomass.native=sum(Biomass[which(is.na(degreeOfEstablishment))]))%>%
	left_join(plots.div[,c("pltID",vars)],by="pltID")
inv.div3=inv.div%>%filter(PLT_CN%in%inv.after.div$PREV_PLT_CN)#390 plots
delta.after=do.call(cbind,Map('-', inv.after.div[,-c(1:3)],inv.div3[,-c(1:3)]))

## t test
dat.pre=rbind(data.frame(invaded=1,PLT_CN=inv.div2$PREV_PLT_CN, inv.div2[,c("SR.native","BA.native","Biomass.native",vars)]),
	data.frame(invaded=0,PLT_CN=inv.pre.div$PLT_CN,inv.pre.div[,c("SR.native","BA.native","Biomass.native",vars)]))
re=c()
for (var.t in c("SR.native","BA.native","Biomass.native",vars)){
	dat=na.omit(dat.pre[,c(var.t,"invaded","PLT_CN")])
	dat[,1]=scale(dat[,1]);colnames(dat)[1]="var.t"	
	a=dat%>%filter(invaded==1);b=dat%>%filter(invaded==0)
	keep=intersect(b$PLT_CN,a$PLT_CN)#native species distinct after Invaded
	dat=dat%>%filter(PLT_CN%in%keep)
	dat$invaded=factor(dat$invaded,levels=c(1,0))
	
	test=t.test(var.t ~ invaded,data=dat,alternative ="two.sided",paired = TRUE, conf.level = 0.95)
	P=test$p.value
	Px=ifelse(P<=0.001,"***",
		ifelse(P>0.001&P<=0.01,"**",
		ifelse(P>0.01&P<=0.05,"*",
		#ifelse(P>0.05&P<0.1,"",
		ifelse(P>0.05,"ns",P))))
	tmp=data.frame(vars=var.t,mean=test$est,se=test$std, P=P,p=ifelse(P<0.001,paste0("<0.001",Px),paste0(round(P,3),Px)))
	re=rbind(re,tmp)
}
#########################################
##### drivers of invasion ###############
#########################################
library(dplyr);library(lme4)
load("data analysis2/plots.summary.rda")
plots=plots.div
plots$fPlot=factor(plots$pltID.no.year)
plots$disturb=factor(plots$disturb,levels=c(0,1))
plots$invaded=factor(plots$invaded,levels=c(0,1))
vars=c("SR","BA","biomass","HillEven","PD","M_Ap","FAD","M_At","Hsp","std")	

#It doesn’t necessarily mess up the fit, but large differences in the scales of parameters often lead to problems (especially with calculating standard deviations of fixed effects); we might as well clean up this problem before we proceed:
#ref:https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
plots[,vars]=scale(plots[,vars])
#Generalized Linear Mixed-Effects Models
formu=as.formula(paste("invaded ~ ",paste(vars,collapse="+"),"+ (1 | ecoregion)+(1|fPlot)",sep=" "))			
M1 <- glmer(formu2,data=plots,family = "binomial")
#M2 <- glmer(invaded ~ SP*BA*biomass*std + disturb + (1 | ecoregion)+(1|fPlot),data=plots,family = "binomial")
save(M1,file="data analysis2/GLMEresults.rda")
dat=as.data.frame(summary(M1)$coef[-1,])
colnames(dat)[c(2,4)]=c("se","p")
dat$var=rownames(dat);dat$var=factor(dat$var,levels=rev(dat$var))
dat$type=ifelse(dat$p>0.05,"insig",ifelse(dat$Estimate>0,"positive","negative"))
ggplot(dat, aes(x=Estimate, y=var,fill=type)) + 	
		#geom_text_repel(aes(x=gradient, y=est,label=est.p,fontface = "italic"),size=2.5)+			
		geom_errorbar(aes(xmin=Estimate-1.96*se, xmax=Estimate+1.96*se,color=type),width=0.3,size=1,show.legend=FALSE,alpha=0.8)+
		scale_fill_manual(values=c("insig"="darkgray", "positive"="blue","negative"="red"))+	
		scale_color_manual(values=c("insig"="darkgray", "positive"="blue","negative"="red"))+	
		geom_point(size=3,shape=21,color="black",show.legend=FALSE) + 
		labs(x="Estimated slope") +
		theme_bw()+
		annotate("text",x=0.5 , y= 1,size=3,
			label = paste0("Plots = ",length(unique(plots$fPlot))),color = "black")+
		geom_vline(xintercept=0,col="darkgray",size=1,linetype="longdash",alpha=0.5)+
		theme(axis.text = element_text(size=12,color='black'),
			axis.title = element_text(size=15,color='black'),
			axis.title.y = element_blank(),
			plot.margin = margin(t = 0.1, r = 0.6, b = 0.1, l = 0.1, "cm"))
#for different year
yr.list=sort(unique(plots$year))
formu2=as.formula(paste("invaded ~ ",paste(vars,collapse="+"),"+ (1 | ecoregion)",sep=" "))	
dat2=c()
for (yr in yr.list){
	plots.yr=plots%>%filter(year==yr)
	M <- glmer(formu2,data=plots.yr,family = "binomial")
	dat.yr=as.data.frame(summary(M)$coef[-1,])
	colnames(dat.yr)[c(2,4)]=c("se","p")
	dat.yr$var=rownames(dat.yr);dat.yr$var=factor(dat.yr$var,levels=rev(dat.yr$var))
	dat.yr$type=ifelse(dat.yr$p>0.05,"insig",ifelse(dat.yr$Estimate>0,"positive","negative"))
	dat.yr$year=yr
	dat2=rbind(dat2,dat.yr)
}
#nplots=plots%>%group_by(year)%>%tally%>%as.data.frame
ggplot(dat2, aes(x=Estimate, y=var,fill=type)) + 	
		#geom_text_repel(aes(x=gradient, y=est,label=est.p,fontface = "italic"),size=2.5)+			
		geom_errorbar(aes(xmin=Estimate-1.96*se, xmax=Estimate+1.96*se,color=type),width=0.3,size=1,show.legend=FALSE,alpha=0.8)+
		scale_fill_manual(values=c("insig"="darkgray", "positive"="blue","negative"="red"))+	
		scale_color_manual(values=c("insig"="darkgray", "positive"="blue","negative"="red"))+	
		geom_point(size=3,shape=21,color="black",show.legend=FALSE) + 
		labs(x="Estimated slope") +
		theme_bw()+		
		geom_vline(xintercept=0,col="darkgray",size=1,linetype="longdash",alpha=0.5)+
		facet_wrap(~year, scales="fixed")+
		theme(axis.text = element_text(size=12,color='black'),
			axis.title = element_text(size=15,color='black'),
			axis.title.y = element_blank(),
			plot.margin = margin(t = 0.1, r = 0.6, b = 0.1, l = 0.1, "cm"),
			strip.text = element_text(colour = 'black', size = rel(1.5)), 
			strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))

##data summary
#invasive richness through time
#no. of plots through time
time.plots=allplots%>%select(pltID,MEASYEAR,invaded)%>%distinct%>%group_by(MEASYEAR)%>%summarize(Invaded=sum(invaded),Uninvaded=(length(pltID)-Invaded))%>%
	tidyr::pivot_longer(cols=Invaded:Uninvaded,names_to = "Type", values_to = "n")

#no.of species richness through time
time.species2.inv=allplots%>%filter(invaded==1)%>%select(Accepted_name,MEASYEAR,degreeOfEstablishment)%>%distinct%>%group_by(MEASYEAR)%>%
	summarize(Invasive=sum(degreeOfEstablishment%in%status[1:2]),#Established=sum(degreeOfEstablishment%in%status[3]),
		Native=sum(is.na(degreeOfEstablishment)))%>%tidyr::pivot_longer(cols=Invasive:Native,names_to = "Type", values_to = "richness")
time.species2.nat=allplots%>%filter(invaded==0)%>%select(Accepted_name,MEASYEAR,degreeOfEstablishment)%>%distinct%>%group_by(MEASYEAR)%>%
	summarize(Type="uninvaded",richness=sum(is.na(degreeOfEstablishment)))
time.species2=rbind(time.species2.inv,time.species2.nat)

#all, inv richness and % of invasive species in each plot
inv.plots=allplots%>%filter(invaded==1)%>%select(MEASYEAR,pltID,Accepted_name,degreeOfEstablishment)%>%distinct%>%group_by(MEASYEAR,pltID)%>%
	summarize(nat.rich=sum(is.na(degreeOfEstablishment)),inv.rich=(sum(degreeOfEstablishment%in%status[1:2])),por=inv.rich/length(degreeOfEstablishment)*100,.groups="keep")
uninv.plots=allplots%>%filter(invaded==0)%>%select(MEASYEAR,pltID,Accepted_name,degreeOfEstablishment)%>%distinct%>%group_by(MEASYEAR,pltID)%>%tally
# mean % of invasive species in each invaded plots
inv.stat=inv.plots%>%group_by(MEASYEAR)%>%summarize(por.mean=mean(por),por.se=sd(por)/sqrt(length(por)))
# total richness  within a plot
all.rich=allplots%>%select(MEASYEAR,pltID,Accepted_name,degreeOfEstablishment,invaded)%>%distinct%>%
	group_by(MEASYEAR,pltID,invaded)%>%tally%>%group_by(MEASYEAR,invaded)%>%summarize(all.mean=mean(n),all.se=sd(n)/sqrt(length(n)),.groups="keep")
all.rich$invaded=factor(all.rich$invaded,levels=c("0","1"))
#mean richness of plots through time
rich.inv=inv.plots%>%group_by(MEASYEAR)%>%summarize(native.mean=mean(nat.rich),inv.mean=mean(inv.rich),native.se=sd(nat.rich)/sqrt(length(nat.rich)),inv.se=sd(inv.rich)/sqrt(length(inv.rich)))
rich.nat=uninv.plots%>%group_by(MEASYEAR)%>%summarize(univ.native.mean=mean(n),univ.native.se=sd(n)/sqrt(length(n)))
native=rich.inv%>%select(MEASYEAR,native.mean,native.se);inv=rich.inv%>%select(MEASYEAR,inv.mean,inv.se)

colnames(rich.nat)=colnames(native)=colnames(inv)=c("MEASYEAR","rich.mean","rich.se")
rich.stat=rbind(data.frame(Type0="Native",Type="Native(Invaded plots)",native),data.frame(Type0="Invasive",Type="Invasive",inv),data.frame(Type0="Native",Type="Native(Uninvaded plots)",rich.nat))

#mean number of occupied plots for each species
occ.inv=allplots%>%select(MEASYEAR,pltID,Accepted_name,degreeOfEstablishment,invaded)%>%distinct%>%
	filter(invaded==1&degreeOfEstablishment%in%status[1:2])%>%group_by(MEASYEAR,Accepted_name)%>%tally()%>%
	left_join(time.plots%>%filter(Type=="Invaded"),by="MEASYEAR")%>%mutate(por=n.x/n.y*100)%>%
	summarize(range.mean=mean(por),range.se=sd(por)/sqrt(length(por)))
occ.nat.invaded=allplots%>%select(MEASYEAR,pltID,Accepted_name,degreeOfEstablishment,invaded)%>%distinct%>%
	filter(invaded==1&is.na(degreeOfEstablishment))%>%group_by(MEASYEAR,Accepted_name)%>%tally()%>%
	left_join(time.plots%>%filter(Type=="Invaded"),by="MEASYEAR")%>%mutate(por=n.x/n.y*100)%>%
	summarize(range.mean=mean(por),range.se=sd(por)/sqrt(length(por)))
occ.nat.uninvaded=allplots%>%select(MEASYEAR,pltID,Accepted_name,degreeOfEstablishment,invaded)%>%distinct%>%
	filter(invaded==0&is.na(degreeOfEstablishment))%>%group_by(MEASYEAR,Accepted_name)%>%tally()%>%
	left_join(time.plots%>%filter(Type=="Uninvaded"),by="MEASYEAR")%>%mutate(por=n.x/n.y*100)%>%
	summarize(range.mean=mean(por),range.se=sd(por)/sqrt(length(por)))	
occ.time=rbind(data.frame(Type="Native(Invaded plots)",occ.nat.invaded),data.frame(Type="Invasive",occ.inv),data.frame(Type="Native(Uninvaded plots)",occ.nat.uninvaded))

## make plots
library(ggpubr)
theme=theme_bw()+
	theme(axis.text = element_text(size=12,color='black'),
		axis.title = element_text(size=15,color='black'),
		legend.position= c(0.65,0.85),
		legend.background = element_rect(fill = NA),
		legend.title=element_blank(),
		legend.text=element_text(face="bold",size=12))
theme.rm=theme(axis.text.x = element_blank(),axis.title.x = element_blank())		
p1=ggplot(rich.stat, aes(x=MEASYEAR, y=rich.mean,fill=Type,shape=Type)) + 	
		geom_errorbar(aes(ymin=rich.mean-1.96*rich.se, ymax=rich.mean+1.96*rich.se,color=Type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
		scale_fill_manual(values=c("#F8766D","#00BA38", "#619CFF"))+	
 		scale_shape_manual(values=c(21,21,24))+			
		#scale_y_continuous(labels = scales::label_number(accuracy = 0.01))+		
		geom_point(size=3,color="black",show.legend=T) + 		
		labs(x="Year",y="Mean species richness within plot") +
		geom_smooth(aes(x = MEASYEAR, y = rich.mean,color=Type),data=rich.stat,method="gam",linewidth=1.5,show.legend=FALSE,se =F,linetype=1)+		
		facet_wrap(~Type0, scales="free_y",ncol = 1)+
		theme+#theme.rm+
		theme(strip.text = element_text(colour = 'black', size = rel(1.5)), 
			strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2),
			legend.position= c(0.65,0.6))
p2=ggplot(occ.time, aes(x=MEASYEAR, y=range.mean,fill=Type,shape=Type)) + 	
		geom_errorbar(aes(ymin=range.mean-1.96*range.se, ymax=range.mean+1.96*range.se,color=Type),width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
		scale_fill_manual(values=c("#F8766D","#00BA38", "#619CFF"))+		
		scale_shape_manual(values=c(21,21,24))+			
		geom_point(size=3,color="black",show.legend=T) + 
		labs(x="Year",y="Mean proportion of \nplots occupied by \nspecies (%)") +
		#geom_smooth(aes(x = MEASYEAR, y =range.mean,color=Type),data=occ.time,method = "gam",size=1.5,show.legend=FALSE,se =F,linetype=1)+		
		theme+theme.rm
p3=ggplot(time.species2,aes(x=MEASYEAR, y=richness,fill=Type))+
	 geom_bar(stat='identity', width=1,colour='black',position='dodge',show.legend=F)+
	scale_fill_manual(values = c("#F8766D","#00BA38", "#619CFF"),labels=c("Invasive","Native(Invaded plots)","Native(Uninvaded plots)"))+
	labs(x="Year",y="Overall species richness")+theme+theme.rm
p4=ggplot(inv.stat, aes(x=MEASYEAR, y=por.mean)) + 	
		geom_errorbar(aes(ymin=por.mean-1.96*por.se, ymax=por.mean+1.96*por.se),width=1,size=0.8,show.legend=FALSE,alpha=0.8)+
		#geom_bar(stat = "identity", width=1,fill="darkgray",color='black') + 
		geom_point(size=3,fill="darkgray",color="black",show.legend=F) + 
		geom_smooth(aes(x = MEASYEAR, y = por.mean),data=inv.stat,method="gam",linewidth=1.5,show.legend=FALSE,se =F,linetype=1,color="black")+			
		labs(x="Year",y="Proportion of \ninvasive species\n in invaded plots (%)") +
		theme	
p5=ggplot(all.rich, aes(x=MEASYEAR, y=all.mean,shape=invaded)) + 	
		geom_errorbar(aes(ymin=all.mean-1.96*all.se, ymax=all.mean+1.96*all.se),color="black",width=0.5,size=1,show.legend=FALSE,alpha=0.8)+
		geom_point(size=3,fill="black",color="black",show.legend=T,alpha=0.8) + 
		scale_shape_manual(values=c(24,21),label=c("Uninvaded","Invaded"))+			
		geom_smooth(aes(x = MEASYEAR, y = all.mean),data=all.rich,method="gam",linewidth=1.5,show.legend=FALSE,se =F,linetype=1,color="black")+			
		labs(x="Year",y="Total species \nrichness of plots") +
		#geom_smooth(aes(x = MEASYEAR, y =range.mean,color=Type),data=occ.time,method = "gam",size=1.5,show.legend=FALSE,se =F,linetype=1)+		
		theme
#ggarrange(p1,p2,p3.1,ncol=1,align = "v",labels="auto")
ggarrange(ggarrange(p1,p5,ncol=1,align = "v",labels="auto",heights=c(1,0.5)),ggarrange(p2,p3,p4,ncol=1,align = "v",labels=c("c","d","e")))

#######################################################################
######## obtain other plot data and invasive species within plots ###########
#######################################################################
	   
get.fia=function(statename){
	require(rFIA)
	require(magrittr)
	require(dplyr)
	state <-  readFIA(paste('data analysis2/FIA20240327/',statename,sep=""),tables=c('PLOT', 'TREE','SURVEY','COND'))
	##  create pltID which can uniquely identify a plot
	## STATECD: State code NUMBER; UNITCD: Survey unit code; COUNTYCD: County code; PLOT: Plot number
	state$TREE$pltID <-  stringr::str_c(state$TREE$UNITCD, state$TREE$STATECD, state$TREE$COUNTYCD, state$TREE$PLOT,state$TREE$INVYR,sep = '_')
	## #### APPLYING FILTERS ####
	## data frames	
	tree_df <- state$TREE %>%
		# records for alive trees only: STATUSCD == 1
		# A tree that has TPA_UNADJ = NA does not contribute to estimates of tree density, abundance, etc. Therefore, these trees should not be included in any diversity index calculations.
		# some NA TPA values  are due to “legacy trees” – those are trees that were measured as part of the older (pre-2000) inventories, which are still monitored to estimate mortality and growth etc. at the population level. Legacy trees have a Y in P2A_GRM_FLG in the TREE table	
		# remove data missing INVYR
		filter(STATUSCD == 1|is.na(STATUSCD))%>% 
		# UNITCD, STATECD,COUNTYCD,PLOT and INVYR will identify a plot record. Most plot records represent measurements. But a few do not. PLOT_STATUS_CD = 1,  will select plot records that correspond to measurements.	
		# ECOSUBCD: classifying the plot by ecoregions
		left_join(select(state$PLOT,MEASYEAR,MEASMON,LAT,LON,ELEV,CN,PLOT_STATUS_CD,SRV_CN,ECOSUBCD), by=c("PLT_CN"="CN")) %>% filter(PLOT_STATUS_CD ==1|is.na(PLOT_STATUS_CD)) %>%
		# select only the data from annual inventories use the SURVEY table: ANN_INVENTORY == Y	
		# the survey is not for a P3 ozone inventory
		left_join(select(state$SURVEY,CN,ANN_INVENTORY,P3_OZONE_IND), by=c("SRV_CN"="CN")) %>% filter((ANN_INVENTORY == "Y"|is.na(ANN_INVENTORY))& (P3_OZONE_IND == "N"|is.na(P3_OZONE_IND))) %>% 
		# Filter out plantations (by set STDORGCD==0); 
		# Filter out STDSZCD 5 (non-stocked plots)
		# COND_STATUS_CD == 1,2 is accessible forest and nonforest land
		# Filter plots with a single condition using a 95% threshold for CONDPROP_UNADJ
		# PHYSCLCD:  coarsely classify plots as to whether they are on xeric (dry), mesic (moist), and hydric (water-logged) soils. 
		left_join(select(state$COND,PLT_CN,STDORGCD,STDSZCD,COND_STATUS_CD,CONDPROP_UNADJ), by="PLT_CN") %>% 
		filter((STDORGCD == 0|is.na(STDORGCD))&(STDSZCD!=5|is.na(STDSZCD))&(CONDPROP_UNADJ >= 0.95|is.na(CONDPROP_UNADJ))&(COND_STATUS_CD%in%c(1,2)|is.na(COND_STATUS_CD))) %>%
		select(pltID,SPCD,MEASYEAR,MEASMON,LAT,LON,ELEV) %>%distinct
	tree_df$statename=statename
	tree_df$Source="FIA"
	return(tree_df)
	rm(state);gc()
}
require(parallel)
no_cores <- detectCores() - 1
mycl <- makePSOCKcluster(no_cores); 
fia=do.call(rbind,parLapply(cl=mycl,state.list,get.fia))
stopCluster(mycl)
sp=read.csv("data analysis2/FIA20240327/REF_SPECIES.csv")
fia=fia%>%left_join(sp[,c("SPCD","SCIENTIFIC_NAME")],by="SPCD")
fia$pltID.no.year = substr(fia$pltID,1,nchar(fia$pltID)-5)

#FIA urban ref:https://www.fs.usda.gov/research/products/dataandtools/tools/urban-datamart
state <-  readFIA('data analysis2/FIA_urban20240327/',tables=c('PLOT', 'TREE','COND'))
statename=data.frame(STATECD=c(6,11,17,18,19,20,23,24,25,27,29,41,44,48,51,53,55),
	statename=c('CA','DC','IL','IN','IA','KS','ME','MD','MA','MN','MO','OR','RI','TX','VA','WA','WI'))
state$TREE$pltID.no.year <-  stringr::str_c(state$TREE$UNITCD, state$TREE$STATECD, state$TREE$COUNTYCD, state$TREE$PLOT,sep = '_')
tree_df <- state$TREE%>%left_join(statename,by="STATECD")%>%select(statename,pltID.no.year,SPCD,STATUSCD,PLT_CN) %>%filter(STATUSCD == 1)%>%
	#Tree without damaging from insects or diseases, etc.
	# records for alive trees only: STATUSCD == 1
	#PLOT_STATUS_CD %in%c(1,2): keep forest and nonforest sample
	left_join(select(state$PLOT,MEAS_YEAR,MEAS_MONTH,LAT,LON,CN,PLOT_STATUS_CD), by=c("PLT_CN"="CN")) %>% filter(PLOT_STATUS_CD %in%c(1,2)) %>%
		# Filter out plantations (by set STDORGCD==0); 
		# COND_STATUS_CD == 1,2 is accessible forest and nonforest land
		# TREATMENT_CD1=0 & DISTURBANCE_CD1=0 to ensure no treatment or disturbance in the plot
		#INVASIVE_STATUS_CD !=3 filter out plots not sampled for invasive plants if they are presented
		left_join(select(state$COND,STDORGCD,COND_STATUS_CD,TREATMENT_CD1,DISTURBANCE_CD1,INVASIVE_STATUS_CD,PLT_CN), by="PLT_CN") %>% 
		filter((is.na(STDORGCD)| STDORGCD== 0)&COND_STATUS_CD%in%c(1,2)&is.na(TREATMENT_CD1)&is.na(DISTURBANCE_CD1)&(is.na(INVASIVE_STATUS_CD)|INVASIVE_STATUS_CD!=3)) %>%
		select(statename,pltID.no.year,SPCD,MEAS_YEAR,MEAS_MONTH,LAT,LON) %>%distinct
tree_df$ELEV=NA
tree_df$Source="FIA_urban"
sp=read.csv("data analysis2/FIA_urban20240327/REF_SPECIES.csv")
tree_df=tree_df%>%left_join(sp[,c("SPCD","SCIENTIFIC_NAME")],by="SPCD")
tree_df$pltID = stringr::str_c(tree_df$pltID.no.year,tree_df$MEAS_YEAR,sep = '_')
tree_df=tree_df%>%select(pltID,SPCD,MEAS_YEAR,MEAS_MONTH,LAT,LON,ELEV,statename,Source,SCIENTIFIC_NAME,pltID.no.year)
colnames(tree_df)=colnames(fia)
fia=rbind(fia,tree_df)
save(fia,file="data analysis2/fia.rda") # elevation in feet
#fia%>%select(MEASYEAR,statename)%>%distinct%>%group_by(MEASYEAR)%>%tally%>%as.data.frame

#vegbank plot
library(data.table);library(dplyr)
sp1=fread("vegbank_export_csv (1)/plot_taxa.csv")
sp1[sp1$currenttaxoninterp_scientificnamenoauthors%in%"","currenttaxoninterp_scientificnamenoauthors"]=sp1[sp1$currenttaxoninterp_scientificnamenoauthors%in%"","currenttaxoninterp_scientificname"]
sp1[sp1$currenttaxoninterp_scientificnamenoauthors%in%"","currenttaxoninterp_scientificnamenoauthors"]=sp1[sp1$currenttaxoninterp_scientificnamenoauthors%in%"","authorplantname_vb"]
sp1=sp1%>%select(observation_id,currenttaxoninterp_scientificnamenoauthors)%>%filter(currenttaxoninterp_scientificnamenoauthors!="")%>%distinct
plot1=fread("vegbank_export_csv (1)/plot_env.csv")%>%select(observation_id,obsstartdate_vb,latitude,longitude,stateprovince_vb,elevation)%>%distinct
vegbank=sp1%>%left_join(plot1,by="observation_id")%>%filter(!is.na(latitude))
vegbank[vegbank$latitude>90,"latitude"]=vegbank[vegbank$latitude>90,"latitude"]-70#a plot in cororado
vegbank[vegbank$longitude>0,"longitude"]=-vegbank[vegbank$longitude>0,"longitude"]# a plot in califonia
vegbank$year=as.numeric(format(vegbank$obsstartdate_vb, "%Y"))
vegbank[vegbank$year==1010,"year"]=2010
vegbank[vegbank$year==3007,"year"]=2007
vegbank[vegbank$year==6006,"year"]=2006
#vegbank%>%select(year,stateprovince_vb)%>%distinct%>%group_by(year)%>%tally%>%as.data.frame
save(vegbank,file="data analysis2/vegbank.rda") # elevation in feet, 1963-2015

load("data analysis2/fia.rda")
load("data analysis2/vegbank.rda")
vegbank$pltID=paste(vegbank$observation_id,vegbank$year,sep="_")
vegbank$Source="vegbank"
vegbank=vegbank%>%select(pltID,observation_id,stateprovince_vb,latitude,longitude,elevation,year,currenttaxoninterp_scientificnamenoauthors,Source)
colnames(vegbank)=c("pltID","pltID.no.year","state","Lat","Lon","elev","year","Spname","Source")
fia=fia%>%select(pltID,pltID.no.year,statename,LAT,LON,ELEV,MEASYEAR,SCIENTIFIC_NAME,Source)
colnames(fia)=colnames(vegbank)
allplots=rbind(fia,vegbank)

spname=unique(allplots$Spname)
library(TNRS)
res=TNRS(taxonomic_names = spname,sources = "wfo")%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1&Overall_score>=0.9)%>%
	dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)
unmat=sort(spname[!spname%in%res$Name_submitted])
write.table(unmat,"unmat.txt")#change Vulpia microstachys var. microstachys to Vulpia microstachys, etc.

# unmat=read.csv("unmat.csv")%>%left_join(read.csv("tnrs_result.csv")[,-1],by=c("Spname2"="Name_submitted"))
# write.table(unmat,"unmat.csv")#manully check, change Spname1 to Name_submitted
res2=read.csv("unmat.csv")%>% dplyr::select(Name_submitted,Taxonomic_status,Accepted_name_id,Accepted_name)%>%rbind(res)
allplots=allplots%>%left_join(res2,by=c("Spname"="Name_submitted"))%>%filter(!is.na(Accepted_name))%>%select(pltID,pltID.no.year,state,Lat,Lon,elev,year,Source,Accepted_name,Accepted_name_id)%>%distinct

#invasive plants
# splist=read.csv("D:/UFpostdoc/data analysis/US-RIIS.csv")%>%filter(phylum%in%"Tracheophyta")%>%select(scientificName,degreeOfEstablishment)
# res=read.csv("data analysis2/tnrs_result_invasive.csv")%>%select(Name_submitted,Accepted_name)%>%left_join(splist,by=c("Name_submitted"="scientificName"))
# write.csv(res,"data analysis2/res.US-RIIS.csv")#manuly check duplicated matches
res=read.csv("data analysis2/res.US-RIIS.csv")%>%select(-Name_submitted)
allplots=allplots%>%left_join(res,by="Accepted_name")
#allplots%>%filter(!is.na(degreeOfEstablishment))%>%select(Accepted_name,degreeOfEstablishment)%>%distinct%>%group_by(degreeOfEstablishment)%>%tally
allplots$elev=allplots$elev*0.3048

status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
inv.plots=allplots%>%filter(degreeOfEstablishment%in%status[1:2])%>%select(pltID)%>%distinct
inv.plots$invaded=1
allplots=allplots%>%left_join(inv.plots,by="pltID")
allplots[is.na(allplots$invaded),"invaded"]=0
save(allplots,file="data analysis2/allplots.rda")
############################################################
######## map plots by intersecting with raster map #########
############################################################
library(dplyr);library(ggpubr)
load("data analysis2/allplots.rda")	
# allplots%>%select(pltID,year,invaded)%>%distinct%>%group_by(year)%>%summarize(Invaded=sum(invaded),Uninvaded=(length(pltID)-Invaded))%>%as.data.frame
# year=sort(unique(allplots$year))
# write.csv(year,"year.csv")
	   
#map all plots and color based on data scource
library(sf);library(raster)
world.rst <- raster("data analysis/wc2.1_10m_bio/wc2.1_10m_elev.tif")
US <- sf::read_sf('data analysis/cb_2018_us_state_20m/cb_2018_us_state_20m.shp')%>%as(., "Spatial")%>%subset(!.@data$NAME%in%c("Hawaii","Alaska","Puerto Rico"))
US.rst=US%>%raster::mask(world.rst,.)%>%raster::crop(extent(US))
US.sf=US.rst%>%stars::st_as_stars()%>%st_as_sf()
plots = allplots%>%dplyr::select(pltID.no.year,Lon,Lat)%>%distinct%>% #filter(Source=="FIA")	
	st_as_sf(coords=c("Lon","Lat"))
st_crs(plots)=st_crs(US.sf)
intersec=st_intersects(US.sf,plots)
domain=raster::as.data.frame(US.rst,xy=TRUE)%>%na.omit
domain$grid.id=1:nrow(domain)
mf=function(i,plots,intersec) if (length(intersec[[i]])>0) data.frame(grid.id=i,pltID.no.year=as.data.frame(plots)[intersec[[i]],-2]) else NA
plots.in.grid=do.call(rbind,lapply(1:length(intersec),mf,plots,intersec))%>%left_join(domain,by="grid.id")%>%na.omit 

year.interval=read.csv("year.csv")
allplots.grid=allplots %>% left_join(plots.in.grid,by="pltID.no.year",relationship ="many-to-many") %>%filter(!is.na(grid.id))%>% left_join(year.interval,by="year")
pltdata=c();pltdata.inv=c()
for (yr in unique(year.interval$interval)){
allplots.grid.yr=allplots.grid%>%filter(interval==yr)%>%dplyr::select(grid.id,Source,interval,invaded)%>%distinct%>%
		right_join(domain,,by="grid.id")%>%as.data.frame
allplots.grid.yr$interval=yr
inv.grid=allplots.grid.yr;
inv.grid[inv.grid$invaded==0&(!is.na(inv.grid$invaded)),"Source"]=NA		
pltdata=rbind(pltdata,allplots.grid.yr)	
pltdata.inv=rbind(pltdata.inv,inv.grid)
}
theme=theme_bw()+
		theme(axis.text = element_text(size=12,color='black'),
			axis.title = element_text(size=15,color='black'),
			strip.text=element_text(size=15,color='black',face="bold"),
		strip.background = element_rect(fill = 'white', colour = 'black', size = rel(2), linetype = 2))	
ggplot(data = pltdata) +
		geom_tile(aes(x = x, y = y,fill = Source)) +
		labs(x="Longitude",y="Latitude")+
		scale_fill_manual(values=c("FIA"="darkred", "vegbank"="blue","FIA_urban"="green"), na.value="lightgray")+
		coord_quickmap() +  theme_bw() + 
		facet_wrap(~interval)+
		theme

################
## data summary#
################
#invasive richness through time
library(dplyr)
#detach(package:raster)
status=c("widespread invasive (category E)","invasive (category D2)","established (category C3)")
load("data analysis2/allplots.rda")
#year.interval=read.csv("year.csv")
#allplots=allplots%>%left_join(year.interval,by="year")%>%filter(Lon>-100&!is.na(year))#%>%filter(Source=="FIA")
allplots$interval=allplots$year
allplots=allplots%>%filter(Lon>-100&interval>=2000&interval<=2021)%>%filter(Source=="FIA")
#no. of plots through time
time.plots=allplots%>%select(pltID,interval,invaded)%>%distinct%>%group_by(interval)%>%summarize(Invaded=sum(invaded),Uninvaded=(length(pltID)-Invaded))%>%
	tidyr::pivot_longer(cols=Invaded:Uninvaded,names_to = "Type", values_to = "n")
p1=ggplot(time.plots, aes(x=interval, y=n,fill=Type)) + 	
		geom_bar(stat="identity",position="stack", color="black", width=1,size=0.25)+
		scale_fill_manual(values = c("#F8766D", "#619CFF"))+
		labs(x="Year",y="Number of plots") +
		theme+theme.rm+theme(legend.position= c(0.8,0.8))
# % of plots contains invasive species
plot.stat=allplots%>%select(pltID,interval,invaded)%>%distinct%>%group_by(interval)%>%summarize(por=sum(invaded)/length(pltID)*100,.groups="keep")
p2=ggplot(plot.stat, aes(x=interval, y=por)) + 	
		geom_bar(stat = "identity", width=1,fill="darkgray",color='black') + 
		labs(x="Year",y="Proportion of \ninvaded plots (%)") +
		theme			
ggarrange(p1,p2,ncol=1,align = "v",labels="auto")

## compare different data sources -----
##No. of plots and species in each data sources
a=allplots%>%select(pltID.no.year,invaded,Source)%>%distinct%>%group_by(Source,invaded)%>%tally
b=allplots%>%select(Accepted_name,degreeOfEstablishment,Source)%>%distinct%>%group_by(Source,degreeOfEstablishment)%>%tally
write.csv(a,"a.csv");write.csv(b,"b.csv")

## statistics for each year
#no of plots in different data scouces
a.inv=allplots%>%filter(invaded==1)%>%select(pltID,Source,year)%>%distinct%>%group_by(year,Source)%>%tally()%>%as.data.frame
a.nat=allplots%>%filter(invaded==0)%>%select(pltID,Source,year)%>%distinct%>%group_by(year,Source)%>%tally()%>%as.data.frame
#no of overall species in different scouces	
b=allplots%>%select(Accepted_name,year,degreeOfEstablishment,Source)%>%distinct%>%group_by(year,Source)%>%
	summarize(Invasive=sum(degreeOfEstablishment%in%status[1:2]),#Established=sum(degreeOfEstablishment%in%status[3]),
		Native=sum(is.na(degreeOfEstablishment)))%>%tidyr::pivot_longer(cols=Invasive:Native,names_to = "Type", values_to = "Richness")
	
p1=ggplot(b,aes(x=year, y=Richness,fill=Source))+
	geom_bar(stat='identity',width=1,colour='black',position='dodge')+
	labs(x="Year",y="Overall species richness")+
	facet_wrap(~Type, scales="free_y",ncol = 1)+
	theme+theme.rm  +theme(legend.position= c(0.8,0.85),strip.text = element_text(colour = 'black', size = rel(1.5)), 
			strip.background = element_rect(fill = 'white', colour = 'darkgray', size = rel(2), linetype = 2))
p2=ggplot(a.inv,aes(x=year,y=n,fill=Source))+
	labs(x="Year",y="Number of \ninvaded plots") +
  geom_bar(stat='identity',width=1,colour='black',position='dodge',show.legend=F)+
  theme+theme() +theme.rm
p3=ggplot(a.nat,aes(x=year,y=n,fill=Source))+
	labs(x="Year",y="Number of \nuninvaded plots") +
  geom_bar(stat='identity',width=1,colour='black',position='dodge',show.legend=F)+
  theme
ggarrange(p1,p2,p3,heights=c(1,0.5,0.5),ncol=1,align="v",labels="auto")
