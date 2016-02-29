#Dot Map Project in R
# an R implementation of https://github.com/unorthodox123/RacialDotMap
# created by: Soon Ju Kim and Justin Joque, University of Michigan

# ***** means possible parallelization

# Initial settings
resolution = function(zoom){
  return(initial.res/(2^zoom))
}

# Read Shapefile, Convert to Coordinates


totalcoordstate2 = function(state, cluster=F) {
  require(sp)
  # state is a spatial polygons df
  withpop = state[state$POP10 != 0, ]
  mapsnpops = cbind(1:nrow(withpop@data),  withpop$POP10)

  if(!(cluster)){
    result = apply(mapsnpops, 1, function(x) spsample(withpop[x[1],], n=x[2], type='random', iter=10))
    pts = sapply(result, coordinates)
  } else {

    require(parallel)
    clus = makeCluster(detectCores()-1)
    clusterExport(clus, c('withpop','mapsnpops'), envir=environment())
    clusterEvalQ(clus, library(sp))

    result = parApply(clus, mapsnpops, 1, function(x) spsample(withpop[x[1],], n=x[2], type='random', iter=10))
    pts = parSapply(clus, result, coordinates)

    stopCluster(clus)
  }

  SpatialPoints(do.call('rbind', pts))
}

# From Lat/Lon to Meters
coordstoMeters2 = function (coords){                                            # ***** (but would need a ridiculous amount of coords for speed gain)
  require(rgdal)                                                                # this projection isn't the same as the shapefile so doesn't work, but the
  coords = sp::SpatialPoints(coords)                                            # approach might be something to consider, as it keeps an SP object.
  proj4string(coords) = CRS("+proj=longlat +ellps=WGS84")
  spTransform(coords, "+proj=utm +ellps=WGS84")
}

# Meters to Pixels
meterstoPixels = function (meters, zoom, origin.shift){
  #calculate resolution, given zoom:                                            # there are pixels and gridtopology classes in sp that might work here,
  res = resolution(zoom)                                                        # but my contextual knowledge is lacking at this point
  data.frame(px = meters$mx / res + origin.shift / res,
             py = meters$my / res + origin.shift / res)
}

# Pixels to Meters
pixelstoMeters = function(px, py, zoom, origin.shift){
	res = resolution(zoom)
	mx = px * res - origin.shift
	my = py * res - origin.shift
	return(c(mx,my))                                                              # not presently used, but is a cbind desired here?
}

# Pixels to Tiles
pixelstoTiles2 = function (pixels, tile.size){
  # pixels is data.frame
  ceiling(pixels/tile.size) - 1
}

# Tile to Quadkey
tilestoQuadkey2 = function(tiles, zoom) {                                       # this function I don't have much context.
  # tiles is data.frame
  #Convert TMS tile coordinates to Quadtree
  quadkey = ""
  tiles[2] = (2 ^ (zoom) - 1) - tiles[2]
  quad = function(zoomlevel){
    digit = 0
    mask = 1 * (2 ^ (zoomlevel - 1)) #what does 1 << (i-1) mean                 # bitwise shift 1 to the left by i-1 amount; in R this is bitwShiftL(a, n)
    if ((bitwAnd(tiles[1], mask) != 0)) {
      digit = digit + 1                                                         # is it desired to have the second if work on the previous?; assuming so.
    }
    if ((bitwAnd(tiles[2], mask) != 0)) {
      digit = digit + 2
    }
    quadkey = paste0(quadkey, digit)
  }
  paste0(sapply(zoom:0, quad), collapse='')                                     # *****  pvec
}

# Stopped here ------------------------------------------------------------


#returns tms value
googleTiles = function(googleTile, level){
  return(c(googleTile[1], 2^level - 1 - googleTile[2]))
}

#quadkey to tile
quadkeytoTiles = function (quadkey){
  tileX = 0
  tileY = 0
  levelofDetail = nchar(quadkey)

  for (i in levelofDetail:1){
    mask = bitwShiftL(1, i - 1)
    char = substr(quadkey, levelofDetail - i + 1, levelofDetail - i + 1)
    if (char == '1') {
      tileX = bitwOr(tileX, mask)
    }
    if (char == '2') {
      tileY = bitwOr(tileY, mask)
    }
    if (char == '3') {
      tileX = bitwOr(tileX, mask)
      tileY = bitwOr(tileY, mask)
    }
  }

  return(c(tileX, tileY, levelofDetail))
}

tilebounds = function(tileX, tileY, levelofDetail,orgiin.shift){
  ll = pixelstoMeters(tileX * tile.size, tileY * tile.size, levelofDetail,origin.shift)
  ur = pixelstoMeters((tileX+1) * tile.size, (tileY+1) * tile.size, levelofDetail,origin.shift)
  return(c(ll[1],ll[2],ur[1],ur[2]))
}

# Draw tiles from quadkey reference
draw.tile =function(quadkey){
  A = 1000
  width = 512
  zoomlevel = nchar(quadkey)
  google_tile = quadkeytoTiles(quadkey)
  tms_tile = googleTiles(google_tile, zoomlevel)
  bounds = tilebounds(tms_tile[1], tms_tile[2], zoomlevel, origin.shift)

  tile_ll = bounds[1] / A
  tile_bb = bounds[2] / A
  tile_rr = bounds[3] / A
  tile_tt = bounds[4] / A

  xscale = width/(tile_rr - tile_ll)
  yscale = width/(tile_tt - tile_bb)
  scale = min(c(xscale,yscale))

  #quad.coord$quadshort = substring(quad.coord$quadkey,1,zoomlevel)
  coords = quad.coord[quad.coord$quadzoom == quadkey,]

  coords$mx = (coords$meters.mx/A - tile_ll) * scale
  coords$my = (coords$meters.my/A - tile_tt) * -scale


  dir.create(paste(zoomlevel, "/", sep = ""), showWarnings = FALSE)
  dir.create(paste(zoomlevel, "/", tms_tile[1], "/", sep=""), showWarnings=FALSE)
  CairoPNG(file = paste(zoomlevel,"/", tms_tile[1], "/", tms_tile[2], ".png", sep=""),
           width=512, height=512, bg = "transparent")
  #plot to the corners
  par(mar=c(0,0,0,0))
  plot(
    coords$mx,
    coords$my,
    pch = 20,
    xlim = c(1, 512),
    ylim = c(512, 1),
    cex = dot_size(zoomlevel),
    col = rgb(0, 0, 0, dot_opac(zoomlevel) / 255),
    axes = FALSE,
    xaxs = "i",
    yaxs = "i",
    bty = "n"
  )
  dev.off()
}

dot_size <- function(zoom) {
  switch(as.character(zoom),
         "4"=0.1,"5"=0.1,"6"=0.1,"7" = 0.1,"8"=0.1,"9"=0.2,"10"=0.2,"11"=0.2,"12"=0.2,"13"=0.2,"14"=0.2)
}

dot_opac <- function(zoom) {
  switch(as.character(zoom),
         "4"=153,"5"=153,"6"=179,"7" = 179,"8"=204,"9"=204,"10"=230,"11"=230,"12"=255,"13"=255,"14"=255)
}
