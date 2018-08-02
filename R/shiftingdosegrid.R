manEMD <- function(structure1, structure2){
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
      fracVector <- directionVector * 0.001 #step length walking in line direction
      steps <- floor(distanceToTraverse /(fracVector[1]^2 + fracVector[2]^2 + fracVector[3]^2)^0.5) #how many steps (how many fracVectors fit inside of distanceToTraverse)
      rise <- rep(NA, steps)
    for (i in 1:steps) {
        alpha <- Apoint + (i-1) * fracVector
        beta <- Apoint + i* fracVector
        rise[i] <- abs(approx3D(janedoe.RTdata$dose, x=(beta[1]), y=(beta[2]), z=(beta[3]), extrapolate = TRUE) - approx3D(janedoe.RTdata$dose, x=(alpha[1]), y=(alpha[2]), z=(alpha[3]), extrapolate = TRUE))
    }
      iterations[i] <- sum(rise, na.rm=TRUE)
    }
    EMDres <- mean(iterations, na.rm = TRUE)
    return(EMDres)
}
 # manEMD(teeth1, teeth2)
    





    # dimnames(newjane)[[1]] <- as.character(as.numeric(dimnames(oldjane)[[1]])+0.01*directionVector[1]) #moves the new dose grid in the dirction + magnitude fo the directionVecotr
    # dimnames(newjane)[[2]] <- as.character(as.numeric(dimnames(oldjane)[[2]])+0.01*directionVector[2])
    # dimnames(newjane)[[3]] <- as.character(as.numeric(dimnames(oldjane)[[3]])+5)
    # approx3D(newjane, x=(Apoint[1]), y=(Apoint[2]), z=(Apoint[3]), extrapolate = TRUE)
    #calculate differece using approx3D, thats the grad. 
    #integral?
   