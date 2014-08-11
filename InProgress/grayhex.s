# File:   	grayhex.s
# Programmers:	Dan Carr & Ru Sun
# Date:		May 3, 1999; Edited Aug 3, 1999
# Purpose:	Draw hexagon pyramid and lop to the peak
#		
# Note:		This could be a bit more auotmatic 

# ________________ Define hexagon pyramid vertices_____________
theta <- 2*pi*(0:5)/6
mat <- cbind(cos(theta),sin(theta),rep(0,6))
peak <- c(0,0,.7)

# find the light direction as 
# the normal to one of the triangles
dir <- rep(0,3)
p1 <- mat[3,]
p2 <- mat[4,]
a <- p2-p1
b <- peak-p2
dir[1] <- a[2]*b[3]-a[3]*b[2]
dir[2] <- a[3]*b[1]-a[1]*b[3]
dir[3] <- a[1]*b[2]-a[2]*b[1]
dir <- dir/sqrt(sum(dir*dir))

# _____________Initialization_________________
norm <- rep(0,3)    # normal vector
cen <- rep(0,3)     # triangle center
intense <- rep(0,7) # light intensity

#______________Calculate intensities for polytope faces____________


for (i in 1:6){
   k <- ifelse(i==6,1,i+1)
   p1 <- mat[i,]
   p2 <- mat[k,]
   a <- p2-p1
   b <- peak-p2
   norm[1] <- a[2]*b[3]-a[3]*b[2]
   norm[2] <- a[3]*b[1]-a[1]*b[3]
   norm[3] <- a[1]*b[2]-a[2]*b[1]
   norm <- norm/sqrt(sum(norm*norm))
   intense[i] <- sum(norm*dir)
}

# set intensity for the top face
intense[7] <- sum(c(0,0,1)*dir)

# __________Scale diffuse intensities and add ambient background________
# based on looking at intensties
# define gray levels
val <- .1+.8*intense
grays <- cbind(val,val,val)
grays <- rbind(c(0,0,0),c(1,1,1),grays,c(.65,.65,.65))
grays <- rgb(grays[,1],grays[,2],grays[,2])

nmat <- rbind(mat,peak)

# _____________Graphics Device and Plot Scaling__________________
#postscript('grayhex.ps',hori=FALSE,width=6.5,height=9,colors=gray)
par(pty='s')
plot(2.5*c(-1,1),2.5*c(-1,1),type='n',xlab='',ylab='')

# ________________Draw triangle faces and lop off peak_____________

nmat <- rbind(mat,peak)
for (i in 1:6){
  k <- ifelse(i==6,1,i+1) 
  subs <- c(i,k,7)
  polygon(nmat[subs,1],nmat[subs,2],col=grays[2+i],border=FALSE)
  #polygon(mat[,1],mat[,2],col="black",density=0)
  # triangle of center hexagon
  polygon(.3*nmat[subs,1],.3*nmat[subs,2],col=grays[9],border=FALSE)
}
#outline center hexagon
polygon(.3*mat[,1],.30*mat[,2],col="black",density=0)
print(mat)
print(nmat)


