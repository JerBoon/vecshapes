
#--------------------------------------------------------------------------------
# Library of various 3D shapes
#--------------------------------------------------------------------------------


#---------------------------------------------------------------------
#' Return a dodecahedron  object, constructed from elementary spacial objects
#'
#' @param centre Centre of dodecahedron 
#' @param radius Is the nominal size from the centre to each vertex. So the radius of the smallest sphere containing the object.
#' @param properties Package-independent object defining additional properties
#' @param bound Include a bounding sphere? Default=TRUE
#' @param align.to Name of axis to align to - i.e. 2 surfaces will be level in that plane. Default = "x"
#'
#' @return A compound object, comprising of elementary triangle objects
#'   in a grouped hierarchical object, with a surrounding bounding sphere
#'   used for spacial indexing.
#' 
#' @export
#' @importFrom vecspace Spc.MakePolygon Spc.Combine Spc.Translate Spc.Rotate Spc.MakeSphere
#'
#' @family constructors
#'
#' @examples
#'   w <- Spc.MakeDodecahedron(c(0,0,0), 10, surface_props)

Spc.MakeDodecahedron <- function (centre, radius, properties=NA, bound=TRUE, align.to="x") {

  if ((typeof(centre) != "double") || (length(centre) != 3)) {
    print("Spc.MakeDodecahedron: centre should be a 3 number vector")
    return(NA)
  }

  if ((typeof(radius) != "double") || (length(radius) != 1)|| (radius <= 0)) {
    print("Spc.MakeDodecahedron: radius should be a positive scalar numeric")
    return(NA)
  }

  #----

  #Definition of vertices, from https://en.wikipedia.org/wiki/Dodecahedron#Cartesian_coordinates

  phi <- (1 +sqrt(5))/2

  verts <- matrix(c(1,1,1,1,1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,
                    0,phi,1/phi,0,phi,-1/phi,0,-phi,1/phi,0,-phi,-1/phi,
                    1/phi,0,phi,1/phi,0,-phi,-1/phi,0,phi,-1/phi,0,-phi,
                    phi,1/phi,0,phi,-1/phi,0,-phi,1/phi,0,-phi,-1/phi,0),
                  ncol=3, byrow=TRUE)

  faces <- matrix(c(9,1,17,2,10,
                    20,8,12,11,7,
                    18,3,11,12,4,
                    10,6,19,5,9,
                    2,17,18,4,14,
                    17,1,13,3,18,
                    6,16,8,20,19,
                    19,20,7,15,5,
                    10,2,14,16,6,
                    13,1,9,5,15,
                    11,3,13,15,7,
                    14,4,12,8,16),
                 ncol=5, byrow=TRUE)

  to.face <- function (v) {
    return(Spc.MakePolygon(verts[v,1] * scale,verts[v,2] * scale,verts[v,3] * scale))
  }

  #The default size is root 3, use in all the calculations of phi above, etc
  #So calculate a scaling factor

  scale <- radius / sqrt(3)

  #----

  r <- list()

  for (i in 1:12)
    r <- append(r, to.face(faces[i,]))

  #Don't use the Spc.Combine bounding algorith, it's naive!
  r <- Spc.Combine(r, properties=properties, bound=FALSE)

  if (bound) {
    sp <- Spc.MakeSphere(c(0,0,0),sqrt(3) * scale)
    sp$objects <- r
    r <- sp
  }

  #---- rotations to plane

  dihangle <- (pi - atan(2)) *180 / pi

  if (align.to %in% c("x","X"))
    r <- Spc.Rotate(r, ,pivot.angle=c(0,0,dihangle/2))

  if (align.to %in% c("y","Y"))
    r <- Spc.Rotate(r, ,pivot.angle=c(dihangle/2,0,0))

  if (align.to %in% c("z","Z"))
    r <- Spc.Rotate(r, ,pivot.angle=c(0,dihangle/2,0))

  r <- Spc.Translate(r, centre)

  return(r)
}





