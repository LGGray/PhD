replace.names <- function(celltype) {
  # Define a list of replacements
  replacements <- list(
    "Age.associated.B.cells" = "Age-associated B cells",
    "all.cells" = "All cells",
    "B.cells" = "B cells",
    "CD16+.NK.cells" = "CD16+ NK cells",
    "Classical.monocytes" = "Classical monocytes",
    "Cycling.T.cells" = "Cycling T cells",
    "DC1" = "DC1", "DC2" = "DC2", "HSC.MPP" = "HSC/MPP",
    "MAIT.cells" = "MAIT cells", "Memory.B.cells" = "Memory B cells",
    "Naive.B.cells" = "Naive B cells", "NK.cells" = "NK cells",
    "Non.classical.monocytes" = "Non-classical monocytes",
    "pDC" = "pDC", "Plasma.cells" = "Plasma cells",
    "Plasmablasts" = "Plasmablasts", "Regulatory.T.cells" = "Regulatory T cells",
    "Tcm.Naive.cytotoxic.T.cells" = "Tcm/Naive cytotoxic T cells",
    "Tcm.Naive.helper.T.cells" = "Tcm/Naive helper T cells",
    "Tem.Effector.helper.T.cells" = "Tem/Effector helper T cells",
    "Tem.Temra.cytotoxic.T.cells" = "Tem/Temra cytotoxic T cells",
    "Tem.Trm.cytotoxic.T.cells" = "Tem/Trm cytotoxic T cells")
  
  # Replace the celltype names with their proper names
  celltype_proper <- unlist(replacements[celltype])
  
  # Return the proper names
  return(celltype_proper)
}