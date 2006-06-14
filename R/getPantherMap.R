"getPantherMap" <-
function(probeID, title, figDir) {
	idx <- c(2,3,4)  #- to plot family, function, process
	load("Z:/rFunction/MousePanther.rd")

	if(missing(figDir)) {
		figDir = "./"
	}
	if( ! file.exists(figDir)) {
		dir.create(figDir, showWarnings = FALSE)
	}

	panther.use <- panther[which(rownames(panther) %in% probeID),]
	name <- colnames(panther)
	for(i in idx) {
		blank <- which(panther.use[,i] == "")
		blank <- c(blank, grep("NOT NAMED", as.vector(panther.use[,i])))
		unclassified <- grep("unclassified", as.vector(panther.use[,i]))

		print(paste("Checking information from: ", name[i], " ...", sep = ""))
		flush.console()
	
		panT <- table(factor(panther.use[-c(blank, unclassified), i]))
		idx2more <- which(panT > 1)

		if(any(length(idx2more)) > 2) {
			filename <- paste(figDir, name[i], "_fdr", title, ".bmp", sep ="")
			bmp(filename = filename, width = 800, height = 600)
			pie(panT[idx2more], main = paste(name[i], " (FDR = ", title, ")", sep = ""))

			totalProbCount <- length(probeID)
			mapProbCount <- dim(panther.use)[1]
			blankProbCount <- length(c(blank, unclassified))
			oneProbCount <- mapProbCount - length(idx2more)

			legend <- paste(paste("Total Probe:  ", totalProbCount, "  Mapped Probe: ", mapProbCount, sep = ""),
							 paste("NoName Probe: ", blankProbCount, "  One per Cat:  ", oneProbCount, sep = ""), sep = "\n")

			mtext(legend, side = 1)

			dev.off()
		}
		else {
			print("... No probes mapped in Panther.")
			flush.console()
		}
	}
}

