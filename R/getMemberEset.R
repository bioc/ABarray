"getMemberEset" <-
function(eset, group, member) {
	pd <- pData(eset)
	grpMember <- pd[, colnames(pd) == group]
	return(eset[, is.element(grpMember, member)])
}

