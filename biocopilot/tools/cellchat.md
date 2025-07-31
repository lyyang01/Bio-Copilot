# CellChat
## load cellchat package in R environment
library(CellChat)
## Create a new CellChat object from a data matrix, Seurat or SingleCellExperiment object
createCellChat(
  object,
  meta = NULL,
  group.by = NULL,
  datatype = c("RNA", "spatial"),
  coordinates = NULL,
  scale.factors = NULL,
  assay = NULL,
  do.sparse = T
)
### the detail description of parameters:
1. object	
a normalized (NOT count) data matrix (genes by cells), Seurat or SingleCellExperiment object
2. meta	
a data frame (rows are cells with rownames) consisting of cell information, which will be used for defining cell groups. If input is a Seurat or SingleCellExperiment object, the meta data in the object will be used
3. group.by	
a char name of the variable in meta data, defining cell groups. If input is a data matrix and group.by is NULL, the input ‘meta' should contain a column named ’labels', If input is a Seurat or SingleCellExperiment object, USER must provide 'group.by' to define the cell groups. e.g, group.by = "ident" for Seurat object
4. datatype	
By default datatype = "RNA"; when running CellChat on spatial imaging data, set datatype = "spatial" and input 'scale.factors'
5. coordinates	
a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
6. scale.factors	
a list containing the scale factors and spot diameter for the full/high/low resolution images.
USER must input this list when datatype = "spatial". scale.factors must contain an element named 'spot.diameter', which is the theoretical spot size; e.g., 10x Visium (spot.size = 65 microns), and another element named 'spot', which is the number of pixels that span the diameter of a theoretical spot size in the original, full-resolution image.
For 10X visium, scale.factors are in the 'scalefactors_json.json'. scale.factors$spot is the 'spot.size.fullres '
7. assay	
Assay to use when the input is a Seurat object. NB: The data in the 'integrated' assay is not suitable for CellChat analysis because it contains negative values.
8. do.sparse	
whether use sparse format
## Add the cell information into meta slot
addMeta(object, meta, meta.name = NULL)
### the detail description of parameters:
1. object	
CellChat object
2. meta	
cell information to be added
3. meta.name	
the name of column to be assigned
## Set the default identity of cells
setIdent(object, ident.use = NULL, levels = NULL, display.warning = TRUE)
### the detail description of parameters:
1. object	
CellChat object
2. ident.use	
the name of the variable in object.meta;
3. levels	
set the levels of factor
4. display.warning	
whether display the warning message
## The ligand-receptor interaction database curated in CellChat tool
CellChatDB.human
## Subset CellChatDB databse by only including interactions of interest
subsetDB(CellChatDB, search = c(), key = "annotation")
### the detail description of parameters:
1. CellChatDB	
CellChatDB databse
2. search	
a character
3. key	
the name of the variable in CellChatDB interaction_input
## Subset the expression data of signaling genes for saving computation cost
subsetData(object, features = NULL)
### the detail description of parameters:
1. object	
CellChat object
2. features	
default = NULL: subset the expression data of signaling genes in CellChatDB.use
## Identify over-expressed signaling genes associated with each cell group
identifyOverExpressedGenes(
  object,
  data.use = NULL,
  group.by = NULL,
  idents.use = NULL,
  invert = FALSE,
  group.dataset = NULL,
  pos.dataset = NULL,
  features.name = "features",
  only.pos = TRUE,
  features = NULL,
  return.object = TRUE,
  thresh.pc = 0,
  thresh.fc = 0,
  thresh.p = 0.05
)
### the detail description of parameters:
1. object	
CellChat object
2. data.use	
a customed data matrix. Default: data.use = NULL and the expression matrix in the slot 'data.signaling' is used
3. group.by	
cell group information; default is 'object@idents'; otherwise it should be one of the column names of the meta slot
4. idents.use	
a subset of cell groups used for analysis
5. invert	
whether invert the idents.use
6. group.dataset	
dataset origin information in a merged CellChat object; set it as one of the column names of meta slot when identifying the highly enriched genes in one dataset for each cell group
7. pos.dataset	
the dataset name used for identifying highly enriched genes in this dataset for each cell group
8. features.name	
a char name used for storing the over-expressed signaling genes in 'object@var.features[[features.name]]'
9. only.pos	
Only return positive markers
10. features	
features used for identifying Over Expressed genes. default use all features
11. return.object	
whether return the object; otherwise return a data frame consisting of over-expressed signaling genes associated with each cell group
12. thresh.pc	
Threshold of the percent of cells expressed in one cluster
13. thresh.fc	
Threshold of Log Fold Change
14. thresh.p	
Threshold of p-values
## Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
identifyOverExpressedInteractions(
  object,
  features.name = "features",
  features = NULL,
  return.object = TRUE
)
### the detail description of parameters:
1. object	
CellChat object
2. features.name	
a char name used for assess the results in 'object@var.features[[features.name]]'
3. features	
a vector of features to use. default use all over-expressed genes in 'object@var.features[[features.name]]'
4. return.object	
whether returning a CellChat object. If FALSE, it will return a data frame containing the over-expressed ligand-receptor pairs
## Compute the communication probability/strength between any interacting cell groups
computeCommunProb(
  object,
  type = c("triMean", "truncatedMean", "thresholdedMean", "median"),
  trim = 0.1,
  LR.use = NULL,
  raw.use = TRUE,
  population.size = FALSE,
  distance.use = TRUE,
  interaction.length = 200,
  scale.distance = 0.01,
  k.min = 10,
  nboot = 100,
  seed.use = 1L,
  Kh = 0.5,
  n = 1
)
### the detail description of parameters:
1. **object**
   - CellChat object
2. **type**
   - Methods for computing the average gene expression per cell group. By default = "triMean", producing fewer but stronger interactions; When setting ‘type = "truncatedMean"', a value should be assigned to ’trim', producing more interactions.
3. **trim**
   - The fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed
4. **LR.use**
   - A subset of ligand-receptor interactions used in inferring communication network
5. **raw.use**
   - Whether use the raw data (i.e., 'object@data.signaling') or the projected data (i.e., 'object@data.project'). Set raw.use = FALSE to use the projected data when analyzing single-cell data with shallow sequencing depth because the projected data could help to reduce the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors.
6. **population.size**
   - Whether consider the proportion of cells in each group across all sequenced cells. Set population.size = FALSE if analyzing sorting-enriched single cells, to remove the potential artifact of population size. Set population.size = TRUE if analyzing unsorted single-cell transcriptomes, with the reason that abundant cell populations tend to send collectively stronger signals than the rare cell populations.
7. **distance.use**
   - Whether use distance constraints to compute communication probability. distance.use = FALSE will only filter out interactions between spatially distant regions, but not add distance constraints.
8. **interaction.length**
   - The maximum interaction/diffusion length of ligands (Unit: microns). This hard threshold is used to filter out the connections between spatially distant regions
9. **scale.distance**
   - A scale or normalization factor for the spatial distances. This values can be 1, 0.1, 0.01, 0.001. We choose this values such that the minimum value of the scaled distances is in [1,2]. When comparing communication across different CellChat objects, the same scale factor should be used. For a single CellChat analysis, different scale factors will not affect the ranking of the signaling based on their interaction strength.
10. **k.min**
    - The minimum number of interacting cell pairs required for defining adjacent cell groups
11. **nboot**
    - Threshold of p-values
12. **seed.use**
    - Set a random seed. By default, set the seed to 1.
13. **Kh**
    - Parameter in Hill function
14. **n**
    - Parameter in Hill function
## Filter cell-cell communication if there are only few number of cells in certain cell groups
filterCommunication(object, min.cells = 10)
### the detail description of parameters:
1. object	
CellChat object
2. min.cells	
the minmum number of cells required in each cell group for cell-cell communication
## Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
computeCommunProbPathway(
  object = NULL,
  net = NULL,
  pairLR.use = NULL,
  thresh = 0.05
)
### the detail description of parameters:
1. object	
CellChat object
2. net	
A list from object@net; If net = NULL, net = object@net
3. pairLR.use	
A dataframe giving the ligand-receptor interactions; If pairLR.use = NULL, pairLR.use = object@LR$LRsig
4. thresh	
threshold of the p-value for determining significant interaction
## Calculate the aggregated network by counting the number of links or summarizing the communication probability
aggregateNet(
  object,
  sources.use = NULL,
  targets.use = NULL,
  signaling = NULL,
  pairLR.use = NULL,
  remove.isolate = TRUE,
  thresh = 0.05,
  return.object = TRUE
)
### the detail description of parameters:
1. object	
CellChat object
2. sources.use, targets.use, signaling, pairLR.use	
Please check the description in function subsetCommunication
3. remove.isolate	
whether removing the isolate cell groups without any interactions when applying subsetCommunication
4. thresh	
threshold of the p-value for determining significant interaction
5. return.object	
whether return an updated CellChat object
## Circle plot of cell-cell communication network
netVisual_circle(
  net,
  color.use = NULL,
  title.name = NULL,
  sources.use = NULL,
  targets.use = NULL,
  idents.use = NULL,
  remove.isolate = FALSE,
  top = 1,
  weight.scale = FALSE,
  vertex.weight = 20,
  vertex.weight.max = NULL,
  vertex.size.max = NULL,
  vertex.label.cex = 1,
  vertex.label.color = "black",
  edge.weight.max = NULL,
  edge.width.max = 8,
  alpha.edge = 0.6,
  label.edge = FALSE,
  edge.label.color = "black",
  edge.label.cex = 0.8,
  edge.curved = 0.2,
  shape = "circle",
  layout = in_circle(),
  margin = 0.2,
  vertex.size = NULL,
  arrow.width = 1,
  arrow.size = 0.2
)
### the detail description of parameters:
1. **net**
   - A weighted matrix representing the connections.
2. **color.use**
   - Colors represent different cell groups.
3. **title.name**
   - The name of the title.
4. **sources.use**
   - A vector giving the index or the name of source cell groups.
5. **targets.use**
   - A vector giving the index or the name of target cell groups.
6. **idents.use**
   - A vector giving the index or the name of cell groups of interest.
7. **remove.isolate**
   - Whether to remove the isolate nodes in the communication network.
8. **top**
   - The fraction of interactions to show.
9. **weight.scale**
   - Whether to scale the weight.
10. **vertex.weight**
    - The weight of vertex: either a scale value or a vector.
11. **vertex.weight.max**
    - The maximum weight of vertex; default = max(vertex.weight).
12. **vertex.size.max**
    - The maximum vertex size for visualization.
13. **vertex.label.cex**
    - The label size of vertex.
14. **vertex.label.color**
    - The color of label for vertex.
15. **edge.weight.max**
    - The maximum weight of edge; default = max(net).
16. **edge.width.max**
    - The maximum edge width for visualization.
17. **alpha.edge**
    - The transparency of edge.
18. **label.edge**
    - Whether or not to show the label of edges.
19. **edge.label.color**
    - The color for single arrow.
20. **edge.label.cex**
    - The size of label for arrows.
21. **edge.curved**
    - Specifies whether to draw curved edges, or not. This can be a logical or a numeric vector or scalar. First the vector is replicated to have the same length as the number of edges in the graph. Then it is interpreted for each edge separately. A numeric value specifies the curvature of the edge; zero curvature means straight edges, negative values means the edge bends clockwise, positive values the opposite. TRUE means curvature 0.5, FALSE means curvature zero.
22. **shape**
    - The shape of the vertex, currently “circle”, “square”, “csquare”, “rectangle”, “crectangle”, “vrectangle”, “pie” (see vertex.shape.pie), ‘sphere’, and “none” are supported, and only by the plot.igraph command. “none” does not draw the vertices at all, although vertex labels are plotted (if given). See shapes for details about vertex shapes and vertex.shape.pie for using pie charts as vertices.
23. **layout**
    - The layout specification. It must be a call to a layout specification function.
24. **margin**
    - The amount of empty space below, over, at the left and right of the plot, it is a numeric vector of length four. Usually values between 0 and 0.5 are meaningful, but negative values are also possible, that will make the plot zoom in to a part of the graph. If it is shorter than four then it is recycled.
25. **vertex.size**
    - Deprecated. Use 'vertex.weight'.
26. **arrow.width**
    - The width of arrows.
27. **arrow.size**
    - The size of arrow.
## Visualize the inferred signaling network of signaling pathways by aggregating all L-R pairs
netVisual_aggregate(
  object,
  signaling,
  signaling.name = NULL,
  color.use = NULL,
  thresh = 0.05,
  vertex.receiver = NULL,
  sources.use = NULL,
  targets.use = NULL,
  idents.use = NULL,
  top = 1,
  remove.isolate = FALSE,
  vertex.weight = 1,
  vertex.weight.max = NULL,
  vertex.size.max = NULL,
  weight.scale = TRUE,
  edge.weight.max = NULL,
  edge.width.max = 8,
  layout = c("circle", "hierarchy", "chord", "spatial"),
  pt.title = 12,
  title.space = 6,
  vertex.label.cex = 0.8,
  alpha.image = 0.15,
  point.size = 1.5,
  group = NULL,
  cell.order = NULL,
  small.gap = 1,
  big.gap = 10,
  scale = FALSE,
  reduce = -1,
  show.legend = FALSE,
  legend.pos.x = 20,
  legend.pos.y = 20,
  ...
)
## Visualization of network using heatmap
netVisual_heatmap(
  object,
  comparison = c(1, 2),
  measure = c("count", "weight"),
  signaling = NULL,
  slot.name = c("netP", "net"),
  color.use = NULL,
  color.heatmap = c("#2166ac", "#b2182b"),
  title.name = NULL,
  width = NULL,
  height = NULL,
  font.size = 8,
  font.size.title = 10,
  cluster.rows = FALSE,
  cluster.cols = FALSE,
  sources.use = NULL,
  targets.use = NULL,
  remove.isolate = FALSE,
  row.show = NULL,
  col.show = NULL
)
## Names of cell states will be displayed in this chord diagram
netVisual_chord_cell(
  object,
  signaling = NULL,
  net = NULL,
  slot.name = "netP",
  color.use = NULL,
  group = NULL,
  cell.order = NULL,
  sources.use = NULL,
  targets.use = NULL,
  lab.cex = 0.8,
  small.gap = 1,
  big.gap = 10,
  annotationTrackHeight = c(0.03),
  remove.isolate = FALSE,
  link.visible = TRUE,
  scale = FALSE,
  directional = 1,
  link.target.prop = TRUE,
  reduce = -1,
  transparency = 0.4,
  link.border = NA,
  title.name = NULL,
  show.legend = FALSE,
  legend.pos.x = 20,
  legend.pos.y = 20,
  nCol = NULL,
  thresh = 0.05,
  ...
)
## Compute and visualize the contribution of each ligand-receptor pair in the overall signaling pathways
netAnalysis_contribution(
  object,
  signaling,
  signaling.name = NULL,
  sources.use = NULL,
  targets.use = NULL,
  width = 0.1,
  vertex.receiver = NULL,
  thresh = 0.05,
  return.data = FALSE,
  x.rotation = 0,
  title = "Contribution of each L-R pair",
  font.size = 10,
  font.size.title = 10
)



