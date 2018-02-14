
# title: "Routing_Project"
# author: "Adithya Ajith"
#date: "7 July 2017"
# output: html_document

library(needs)
needs(rjson , RCurl , googleway , dplyr, osrm ,leaflet)


src = c(30.538346,-96.691991) #(lat,lon)
dst = c(30.687838,-96.372954) #(lat,lon)

##Routing with OSRM

req = paste("http://router.project-osrm.org/", "route/v1/", "driving", "/", src[2], ",", src[1], ";", dst[2], ",", dst[1], 
            "?alternatives=false&geometries=polyline&steps=true&overview=full&annotations=true", sep = "")
resRaw = RCurl::getURL(utils::URLencode(req), useragent = "R-User")
vres = jsonlite::validate(resRaw)[1]
if (!vres)
  resRaw = gsub(pattern = "[\\]", replacement = "tilapia", x = resRaw)
route_json = jsonlite::fromJSON(resRaw)
if (!vres)
  route_json$routes$geometry = gsub(pattern = "tilapia", replacement = "\\\\", x = route_json$routes$geometry)
rm(vres,resRaw,req,src,dst)

##Google elevation service
key_g_elevation = "key"

elevation_json = googleway::google_elevation(polyline = route_json$routes$geometry, key = key_g_elevation , location_type = "path" , samples = 512)
ele_dump = cbind(elevation_json$results$location , elevation_json$results[c("elevation","resolution")] )
colnames(ele_dump)[2] = "lon"
rm(elevation_json,key_g_elevation)
plot(x= ele_dump$lon , y= ele_dump$lat)
plot.ts(ele_dump$elevation)
#write.csv( file="route.csv", x= ele_dump , row.names = FALSE )

#tabulating intersection on the route
intersections=NULL
for (i in 1:(nrow(route_json$routes$legs[[1]][1,1][[1]])-1))
intersections=rbind(intersections,route_json$routes$legs[[1]][1,1][[1]]["intersections"][i,][[1]])
intersections= intersections[2:nrow(intersections),]
coor=as.data.frame(matrix(unlist(intersections["location"]), nrow = nrow(intersections["location"]), byrow = T))
colnames(coor)=c("lon","lat")
intersections=cbind(intersections[c("out","entry","bearings","in")] , coor[c("lat", "lon")] )
#write.csv( file="intersections.csv", x= coor[c("lat","lng")] , row.names = FALSE )
rm(coor, i)

#nodes
nodes=data.frame(node=route_json$routes$legs[[1]]$annotation$nodes[[1]],dist_prv_node=c(0,route_json$routes$legs[[1]]$annotation$distance[[1]]) , intersection = NA  , way_id =NA)
nodes$dist_src=0
for(i in 2:nrow(nodes))
  nodes$dist_src[i]=nodes$dist_prv_node[i]+nodes$dist_src[i-1]
nodes=cbind(nodes,googleway::decode_pl(route_json$routes$geometry))
rm(i,route_json)
for( i in  1: nrow(intersections)){ #identifying which nodes are intersectinos too.
  dist = (nodes$lat - intersections$lat[i])^2 + (nodes$lon - intersections$lon[i])^2
  nodes$intersection[which.min(dist)]= i
  }
rm(dist,i)
#write.csv( file="nodes.csv", x= nodes , row.names = FALSE )


## Turning Angle Calculation

anglefun <- function(xx,yy,bearing=TRUE,as.deg=FALSE){
  ## calculates the compass bearing of the line between two points
  ## xx and yy are the differences in x and y coordinates between two points
  ## Options:
  ## bearing = FALSE returns +/- pi instead of 0:2*pi
  ## as.deg = TRUE returns degrees instead of radians
  c = 1
  if (as.deg){
    c = 180/pi
  }
  
  b<-sign(xx)
  b[b==0]<-1  #corrects for the fact that sign(0) == 0
  tempangle = b*(yy<0)*pi+atan(xx/yy)
  if(bearing){
    #return a compass bearing 0 to 2pi
    #if bearing==FALSE then a heading (+/- pi) is returned
    tempangle[tempangle<0]<-tempangle[tempangle<0]+2*pi
  }
  return(tempangle*c)
}

bearing.ta <- function(loc1,loc2,loc3,as.deg=T){
  ## calculates the bearing and length of the two lines formed by three points
  ## the turning angle from the first bearing to the second bearing is also calculated
  ## locations are assumed to be in (X,Y) format.
  ## Options:
  ## as.deg = TRUE returns degrees instead of radians
  if (length(loc1) != 2 | length(loc2) != 2 | length(loc3) !=2){
    print("Locations must consist of either three vectors, length == 2,
or three two-column dataframes")
    return(NaN)
  }
  c = 1
  if (as.deg){
    c = 180/pi
  }
  
  locdiff1<-loc2-loc1
  locdiff2<-loc3-loc2
  bearing1<-anglefun(locdiff1[1],locdiff1[2],bearing=F)
  bearing2<-anglefun(locdiff2[1],locdiff2[2],bearing=F)

  if(is.data.frame(locdiff1)){
    dist1<-sqrt(rowSums(locdiff1^2))
    dist2<-sqrt(rowSums(locdiff2^2))
  }else{
    dist1<-sqrt(sum(locdiff1^2))
    dist2<-sqrt(sum(locdiff2^2))
  }
  
  ta=(bearing2-bearing1)
  
  ta[ta < -pi] = ta[ta < -pi] + 2*pi
  ta[ta > pi] = ta[ta > pi] - 2*pi
  return(list(bearing1=unlist(bearing1*c),bearing2=unlist(bearing2*c),
  ta=unlist(ta*c),dist1=unlist(dist1),dist2=unlist(dist2)))
}

#Turing Angles for equidistant points
theta=numeric()
for(i in 1:(nrow(ele_dump) - 2)){
  a=bearing.ta(ele_dump[i,c("lon","lat")], ele_dump[i+1, c("lon","lat")] , ele_dump[i+2 , c("lon","lat")])  #should pass (longitude,lattitude)
  theta[i]=a$ta
  }
#plot.ts(theta , type="S")
ele_dump$theta = c(NA,theta,NA)

plotrix::twoord.plot(ele_dump$lon , ele_dump$lat ,ele_dump$lon, ele_dump$theta , xlab="Longitude" , ylab="Latitude" , rylab="Angle in Degrees" ,lcol=4 , type=c("l" , "l") , main = "Turning Angle at each Equidistant Point" , lpch = 23)

#Turing angle for nodes
theta=numeric()
for(i in 1:(nrow(nodes) - 2)){
  a=bearing.ta(nodes[i,c("lon","lat")], nodes[i+1, c("lon","lat")] , nodes[i+2 , c("lon","lat")])  #should pass (longitude,lattitude)
  theta[i]=a$ta
  }
plot.ts(theta , type="h")
nodes$theta = c(NA,theta,NA)

plotrix::twoord.plot(nodes$lon , nodes$lat ,nodes$lon, nodes$theta , xlab="Longitude" , ylab="Latitude" , rylab="Angle in Degrees" ,lcol=4 , type=c("l" , "l") , main = "Turning Angle at each Node" , lpch = 23)

rm(anglefun, bearing.ta, a, i, theta)


## fetching OSM Data

ways=data.frame((matrix(ncol = 3, nrow = 0)) )
colnames(ways)=c("way_id","prime_nodes", "offroad")

start.time <- Sys.time()
for( i in 1:nrow(nodes)){
print(i)
req = paste0("http://api.openstreetmap.org/api/0.6/node/",nodes$node[i],"/ways")
resRaw = RCurl::getURL(utils::URLencode(req), useragent = "R-User")
xml_data <- XML::xmlToList(resRaw)
for(j in 1:(length(xml_data)-1)){
a=xml_data[[j]]

if(!a$.attrs["id"] %in% ways$way_id){
   ways[nrow(ways)+1, "way_id"] = a$.attrs["id"]
   ways[nrow(ways),"prime_nodes"]= list(nodes$node[i])
   ways[nrow(ways),"offroad"] = TRUE
   tag = names(a)=="tag"
   for(n in 1:length(tag))
     if(tag[n])
      ways[nrow(ways), a[[n]]["k"]] = a[[n]]["v"]
}
else{
  ways$offroad[ways$way_id %in%  a$.attrs["id"]] = FALSE
  ways$prime_nodes[ways$way_id %in%  a$.attrs["id"]]  =list(c(unlist(ways$prime_nodes[ ways$way_id %in%  a$.attrs["id"]]), nodes$node[i]))
}
}
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


for(q in 1:nrow(ways))  #adding the way information to the nodes along the route
  if(ways$offroad[q] == FALSE)
    nodes$ways[nodes$node %in% unlist(ways$prime_nodes[q])] = q
    
    
rm(end.time,start.time , i , j, n ,resRaw , a, xml_data, q, tag , req)

#write.csv( file="ways.csv", x= ways[, ! colnames(ways) %in% c("prime_nodes")] , row.names = FALSE )

###  Mapping with leaflet

#Marking the nodes and intersection on OSM
nodes_map = leaflet(nodes) %>% 
  addTiles() %>% 
  fitBounds(min(nodes$lon), min(nodes$lat), max(nodes$lon), max(nodes$lat)) %>%
  addCircles(~lon, ~lat, label=paste(paste0("# =" , 1:nrow(nodes)),paste0("Turning Angle = ", as.character(round(nodes$theta,4))),paste0("Highyway_type = ", as.character(ways$highway[nodes$ways])), paste0("#Lanes = ",as.character(ways$lanes[nodes$ways])),paste0("MaxSpeed = ", as.character(ways$maxspeed[nodes$ways])),paste0("Oneway = ",as.character(ways$oneway[nodes$ways])),sep=" || "), weight = 5, radius=5,color=ifelse(is.na(nodes$intersection),"blue","red"), stroke = T, fillOpacity = 0.8)%>% 
  addLegend(position = 'bottomright', colors = c("red","blue"), labels = c("Intersection","Non-Intersection"), opacity = 0.4,
                      title = 'Labels') %>% 
  addProviderTiles(providers$OpenStreetMap)
nodes_map
htmlwidgets::saveWidget(nodes_map, file="nodes.html")

#Marking the equidistant points on OSM
ele_map = leaflet(ele_dump) %>% 
  addTiles() %>% 
  fitBounds(min(ele_dump$lon), min(ele_dump$lat), max(ele_dump$lon), max(ele_dump$lat)) %>%
  addCircles(~lon, ~lat, label=paste0("Turning_Angle = ",as.character(round(ele_dump$theta,4))) , weight = 5, radius=5,color="blue", stroke = TRUE, fillOpacity = 0.8) %>% 
  addProviderTiles(providers$OpenStreetMap)
ele_map
htmlwidgets::saveWidget(ele_map, file="equi_map.html")