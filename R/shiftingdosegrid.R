manEMD <- function(structure1, structure2, doseGrid){ 
    structure1weight <- matrix(c(rep(1, times = length(structure1)/3)))
    structure2weight <- matrix(c(rep(1, times = length(structure2)/3)))
    where.on.structure1 <- attributes(emdw(structure1, structure1weight, structure2, structure2weight, flows = TRUE))$flows[[1]] + 1;
    where.on.structure2 <- attributes(emdw(structure1, structure1weight, structure2, structure2weight, flows = TRUE))$flows[[2]] + 1;
    iterations <- rep(NA, length(where.on.structure1))
    for (i in 1:length(where.on.structure1)) {
      Apoint <- structure1[where.on.structure1[i],]
      Bpoint <- structure2[where.on.structure2[i],]
      distanceToTraverse <- ((Bpoint[1]-Apoint[1])^2 + (Bpoint[2]-Apoint[2])^2 + (Bpoint[3] - Apoint[3])^2)^.5 #total distance between starting and ending points
      distTrav <- c((Bpoint[1]-Apoint[1]), (Bpoint[2]-Apoint[2]), (Bpoint[3] - Apoint[3])) #FIX DIVIDE BY ZERO ERROR FOR STEPS
      directionVector <- c(Bpoint[1]-Apoint[1], Bpoint[2]-Apoint[2], Bpoint[3]-Apoint[3]) / distanceToTraverse #unit length vector pointing in direction of line
      fracVector <- directionVector * 0.1 #step length walking in line direction
      steps <- floor(distanceToTraverse /(fracVector[1]^2 + fracVector[2]^2 + fracVector[3]^2)^0.5) #how many steps (how many fracVectors fit inside of distanceToTraverse)
      if(distanceToTraverse !=  0){
        rise <- rep(NA, steps)
      }else{
        iterations[i] <- 0
        next
      }
    for (j in 1:steps) {
        alpha <- Apoint + (j-0.5) * fracVector
        beta <- Apoint + (j+0.5 )* fracVector
        rise[j] <- abs(approx3D(doseGrid, x=(beta[1]), y=(beta[2]), z=(beta[3]), extrapolate = TRUE) - approx3D(doseGrid, x=(alpha[1]), y=(alpha[2]), z=(alpha[3]), extrapolate = TRUE))
    }
      iterations[i] <- sum(rise, na.rm=TRUE)
    }
    EMDres <- mean(iterations, na.rm = TRUE)
    return(EMDres)
}

 






   